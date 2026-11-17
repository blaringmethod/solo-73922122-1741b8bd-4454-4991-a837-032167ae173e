# Paired FASTQ Re-pairing and Adapter-Trimming Service

测序湖入湖工具：接收顺序被打乱的双端未压缩 FASTQ（R1/R2），按读段号重新
配对，按 R1 顺序交付两份文件，原文逐字节保留。配对成功的任务还可以再按
1～4 条接头序列裁剪 3′ 端接头污染，结果以 ZIP 交付（内含保留/丢弃对数）。
FastAPI + SQLite 实现，单文件上限 256 MiB，进程内存上限 256 MiB，可处理
总量超过内存的输入。

## 配对规则

- 读段号取头行首个空白前的字段，去掉开头 `@`，再去掉末尾的 `/1` 或 `/2`
  （例如 `@read-7/1 sample=x` → `read-7`）。
- 输出按 **R1 的原始顺序**排列；R2 按相同读段号对齐到同一位置。
- 四行记录原文（含头行描述、换行风格）原样写入结果，不做任何改写。

## 失败即整次失败，不开放部分结果

任一情况发生，任务整体标记为 `failed`/`cancelled`，已写入的暂存全部删除，
下载接口返回 409：

- 同一文件内读段号重复（含 `/1`、`/2` 归一化后重复）
- 孤端（某读段号只出现在一侧）
- 格式错误：非四行式、头行不以 `@` 开头、第三行不以 `+` 开头、空文件
- 序列与质量字符串长度不等

## 内存与大输入

- 上传经 `SpooledTemporaryFile` 落盘（超过 1 MiB 即写临时文件），分块复制，
  单文件流式校验 256 MiB 上限，超限返回 413。
- 解析按 512 条/批写入 SQLite（WAL、小 page cache、临时表走文件
  `temp_store=FILE`、关闭 mmap），任意时刻内存里只有一个批次。
- 配对发布使用 `JOIN ... ORDER BY a.ord` 流式游标，分批 `fetchmany` 写出。
- 因此输入总量可以远超 256 MiB；自测以约 393 MiB 总量（两份各 196.5 MiB、
  各 200 万条）验证，峰值 RSS 约 19 MiB（见下文自测结果）。

## 任务与取消语义

- `POST /tasks` 上传后立即返回 `202`，后台守护线程处理。
- `GET /tasks/{id}` 查询状态、阶段（upload/parse/validate/publish/done）与
  0–1 进度。
- `POST /tasks/{id}/cancel`：置取消标记后**立即返回**，不等待配对结束；
  工作线程在下一个批次边界中止。取消与完成之间用条件 UPDATE 仲裁：
  - 取消先到 → 任务最终 `cancelled`，结果永不发布，暂存目录删除；
  - 完成先到 → 取消请求得到 `409`，明确告知任务已完成、不可取消，结果保留。
- 取消后可立刻创建新任务；任务目录以 UUID 隔离，清理只删除本任务目录，
  不会误删新任务文件。
- 服务重启（启动钩子）：所有非终态任务标记 `failed`（原因
  `interrupted by service restart`）并清理暂存与记录；`succeeded` 任务的
  结果文件保留，仍可下载。未完成注册的裁剪目录（进程在打包中途退出）
  也会被清理；已登记为成功的裁剪 ZIP 保留。

## API

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/tasks` | multipart 字段 `r1`、`r2`，返回任务对象（202） |
| GET  | `/tasks/{id}` | 状态/进度/消息/下载链接 |
| POST | `/tasks/{id}/cancel` | 请求取消（立即返回） |
| GET  | `/tasks/{id}/download/r1` | 下载配对后的 R1 |
| GET  | `/tasks/{id}/download/r2` | 下载配对后的 R2 |
| POST | `/tasks/{id}/trims` | 对成功任务裁剪接头，返回 trim 元数据（201） |
| GET  | `/tasks/{id}/trims/{trim_id}` | 查询某次裁剪的元数据与计数 |
| GET  | `/tasks/{id}/trims/{trim_id}/download` | 下载裁剪 ZIP |
| GET  | `/health` | 健康检查 |

启动后打开交互式文档：<http://localhost:8000/docs>，可直接在页面上演示
上传、取消、下载、裁剪全流程。

## 接头裁剪

```
POST /tasks/{task_id}/trims
{"adapters": ["AGATCGGAAGAGCACACGTCTGAACTCCAGTC"]}
```

路径中的任务号即成功任务号；请求体可选带 `task_id`，若带则必须与路径一致
（不一致返回 422）。

- 仅接受**配对成功**（`succeeded`）的任务；其他状态返回 409。
- 接头：1～4 条，每条 16～32 位，仅 `A/C/G/T`（大小写均可，按大写处理）。
  不合法返回 422。
- 读长：每端 40～1000 碱基。源文件中出现范围外的读段时整次 409（不会给出
  半成品）；源文件超过 200,000 对时整次 413。
- 命中规则：读段 3′ 后缀与某条接头**前缀**对齐，重叠 ≥ 8 位且最多 1 处
  错配；多条接头/多个重叠同时命中时取**最长**重叠，其它较短命中忽略。
- 序列与质量字符串在**同一位置**同步截短（头行、`+` 行原样保留，换行风格
  保持）；若任一端裁剪后不足 20 碱基，则**整对丢弃**。
- 其余读对严格按 **R1 原始顺序**写入 ZIP 的 `R1.fastq` / `R2.fastq`，另附
  `summary.json`，`retained_pairs` / `dropped_pairs` / `total_pairs` 即保留与
  丢弃对数。响应体与 `/tasks/{id}/trims/{trim_id}` 也返回相同计数。
- 读取源文件时会复核：双端读号一致、四行结构与 `@`/`+` 起始、序列与质量
  长度一致。文件缺件、尾部截断、读号错配等任一问题都会使**整次请求失败**
  （409/500），不会产生或返回部分 ZIP。
- 原任务结果文件只读，不会被修改，也不会被本次请求删除。
- 每次请求使用独立的 `tasks/<id>/trims/<trim_uuid>/` 目录：同一任务的不同
  接头请求并行时，一次失败或中断只清理自己的目录，不会串换、删除、覆盖
  另一请求的结果。ZIP 先写临时名，完整落盘后才改名，失败即整体删除。
- 裁剪按流式读取/写出（约 1 MiB 批次）；测试覆盖 20 万对、每端 400 碱基
  （约 330 MiB）的源文件，子进程峰值 RSS 约 17 MiB。

示例（curl）：

```bash
TASK=<task id>
curl -s -X POST "http://localhost:8000/tasks/$TASK/trims" \
  -H 'Content-Type: application/json' \
  -d '{"adapters":["AGATCGGAAGAGCACACGTCTGAACTCCAGTC"]}'
# 用返回的 download_url 下载：
curl -o trimmed.zip "http://localhost:8000/tasks/$TASK/trims/<trim_id>/download"
```

## 启动后的两段流程

1. `POST /tasks` 上传乱序双端 FASTQ，轮询 `GET /tasks/{id}` 到 `succeeded`。
2. `POST /tasks/{id}/trims` 提交接头列表，拿到 `download_url` 后下载 ZIP。

## 运行（Docker Compose，内存上限 256 MiB）

```bash
docker compose up --build
# 打开 http://localhost:8000/docs
```

`docker-compose.yml` 中 `mem_limit: 256m` 即交付要求的内存上限；数据保存在
命名卷 `fastq-data`（容器内 `/data`，含 SQLite 与各任务目录），重启后成功
结果仍可下载。

## 本地直接运行（开发/自测）

```bash
python3 -m venv --without-pip .venv
.venv/bin/python tools/get-pip.py
.venv/bin/pip install -r requirements.txt
FASTQ_DATA_DIR=./data .venv/bin/uvicorn app.main:app --port 8000
```

## 演示脚本

```bash
./scripts/demo.sh                 # 需要一个运行在 localhost:8000 的服务
BASE_URL=http://localhost:8000 ./scripts/demo.sh
```

脚本会生成一对乱序 FASTQ（每端 40 碱基），演示上传 → 轮询进度 → 下载校验，
再演示一次带接头的裁剪（含一对被整对丢弃）、一次取消和一次错误输入。

## 自测

```bash
# 短测试（功能/取消/重启语义 + 取消与完成的仲裁 + 接头裁剪），共 112 项检查
bash scripts/run_tests.sh
# 或分别运行
.venv/bin/python -I tests/test_api.py          # 30 项：配对、整次失败、取消、重启恢复
.venv/bin/python -I tests/test_cancel_race.py  # 9 项：取消/完成条件 UPDATE 仲裁
.venv/bin/python -I tests/test_trim.py         # 73 项：裁剪算法、API、校验、并发隔离、限额

# 配对内存压力测试：两份各 196.5 MiB（合计约 393 MiB，超过 256 MiB 内存上限），
# 子进程实测峰值 RSS 约 19 MiB，并逐条校验配对顺序与原文
.venv/bin/python -I tests/test_memory.py 2000000

# 裁剪内存压力测试：20 万对、每端 400 碱基（约 330 MiB），峰值 RSS 约 17 MiB，
# 并流式校验保留顺序、丢弃规则与 summary 计数
.venv/bin/python -I tests/test_trim_memory.py
```

## 目录结构

```
app/
  config.py    # 上限、路径、批大小、裁剪参数等配置（环境变量可覆盖）
  db.py        # SQLite schema、迁移、任务/记录/裁剪访问
  fastq.py     # 流式解析、校验、JOIN 配对发布
  trim.py      # 接头匹配、流式裁剪、ZIP 打包与每请求目录隔离
  worker.py    # 后台线程、取消仲裁、启动恢复、任务目录清理
  main.py      # FastAPI 路由、有界流式上传、裁剪 API
tests/         # 功能自测、裁剪自测与内存压力测试
scripts/demo.sh, scripts/run_tests.sh
Dockerfile, docker-compose.yml, requirements.txt
```

## 遗留项 / 边界说明

- 仅接受未压缩 FASTQ；`.gz` 等压缩格式不在范围内。
- 裁剪请求是同步完成的（HTTP 请求等待 ZIP 生成）；20 万对实测约 25 秒。
  响应中不提供进度轮询，重试会生成一份新的独立 ZIP（不同 `trim_id`）。
- 成功任务与裁剪结果都没有内置 TTL 清理，会一直保留在数据卷中，需外部
  运维定期清理。
- 服务为单实例 + 单 SQLite 文件设计；多实例横向扩展需要换用共享存储/队列。
- 开发沙箱内未安装 Docker，`docker compose up` 未实跑；compose/Dockerfile 仅做
  了 YAML 与结构校验，应用本身已用等价的本地 uvicorn 命令完整验证（含
  `/docs`、上传、取消、下载、重启恢复、裁剪）。在装有 Docker 的环境执行
  `docker compose up --build` 即可，`mem_limit: 256m` 即内存上限。

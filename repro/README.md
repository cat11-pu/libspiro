# 退化输入加固：复现、修复与边界扫描

本目录针对 libspiro 在**退化控制点、重复点、超大坐标**下产生 NaN/Inf 或
崩溃的问题，提供独立于 autotools 的复现与验证工具（不改动原有构建）。

## 1. 问题根因（修复前）

全部问题最终都把非有限值送进了用户的 `bezctx` 回调（moveto/lineto/
quadto/curveto），而公开 API 仍返回成功 `1`。

1. **零缩放比例导致除零（重复点/退化开放轮廓）**
   `setup_path0()` 计算包围盒时，开放轮廓的扫描上界是 `n_seg`，
   而开放轮廓 `n_seg == n-1`，因此**永远漏掉最后一个端点**。当轮廓的全部
   跨度都落在末端点（例如末端是 `v`/`}`，或多个点重合）时，跨度为 0，
   `dm[0]`（scale）变成 0，随后
   `r[i].x = (src[i].x - dm[1]) / dm[0]` 除零得到 NaN，并一路传播到
   `lineto(-nan,-nan)`。

2. **坐标跨度溢出（超大坐标，如 ~1e308）**
   跨度 `xmax-xmin` 超过 double 上限时 `dm[0] = inf`，归一化坐标变成
   `inf/inf = NaN`，求解器在第 1 次迭代即失败；部分场景在曲线输出阶段
   `scale*cos(...)`、`x*dm[0]+dm[1]` 溢出为 ±inf 后仍被发出。

3. **有限性检查默认关闭、输出端无闸门**
   输入/弦长的 `IS_FINITE` 检查被 `#ifdef CHECK_INPUT_FINITENESS` 包裹，
   而该宏默认（`./configure` 不带 `--enable-test_inputs`）不定义；
   `spiro_to_bpath0()` 与 `spiro_seg_to_bpath1()` 对求解结果和最终坐标
   都不做有限性校验，NaN/Inf 直接泄漏到回调。

## 2. 修复方案（全部在 `spiro.c`，合法轮廓逐位不变）

- **输入有限性检查无条件生效**：`setup_path0()` 开头拒绝任何 NaN/Inf 的
  `x/y`；弦长 `hypot` 有限性检查也无条件化（不再依赖编译开关）。
- **零跨度安全钳制**：完全保留历史包围盒算法；仅当算出的 `dm[0]==0`
  时，补扫包含末端点的范围重新计算；若跨度确实为 0（点全部重合），把
  scale 钳制为安全正值 `1.0`，杜绝除零。非退化轮廓走原路径，结果不变。
- **溢出显式报错**：重算后若 `dm[0]` 仍非有限（跨度溢出 double），
  `setup_path0()` 释放内存并返回 `NULL`，上层 API 返回 `0`（明确失败）。
- **输出端有限性闸门**：
  - `spiro_seg_to_bpath1()` 由 `void` 改为返回 `int`（1 成功 / 0 出现
    非有限值），在每个除法（`seg_ch/ch`，含 `ch==0` 守卫）、细分中点和
    每个发往 bezctx 的坐标前做 `IS_FINITE` 校验，递归调用传播失败。
  - `spiro_to_bpath0()` 在构建前校验 `dm[]`、所有段的 `x/y/seg_ch/
    seg_th/ks[0..3]` 以及 moveto 坐标；任一段返回 0 即整体返回 0。

契约：任何退化输入现在要么**返回明确错误（API 返回 0，不产生回调）**，
要么被**安全钳制后输出有限坐标（API 返回 1）**；保证零 NaN/Inf、零崩溃。

## 3. 复现与验证

```sh
sh repro/build.sh        # 构建 repro_default（默认配置）与 repro_checked
./repro/repro_default    # 固定用例：修复后 non_finite_outputs 全为 0

# 边界扫描（含 ASan+UBSan），456 个结构化用例
cc -g -O1 -fsanitize=address,undefined -fno-sanitize-recover=all \
   -I. -Irepro spiro.c bezctx.c spiroentrypoints.c repro/scan.c -lm \
   -o repro/scan_asan
./repro/scan_asan

# 随机 fuzz（含 ASan+UBSan），5 万~10 万次
cc -g -O1 -fsanitize=address,undefined -fno-sanitize-recover=all \
   -I. -Irepro spiro.c bezctx.c spiroentrypoints.c repro/fuzz.c -lm \
   -o repro/fuzz_asan
./repro/fuzz_asan 50000

# 合法轮廓逐位回归（与 git HEAD 原版对比，需先用原版生成 /tmp/orig.txt）
cc -O2 -I. -Irepro spiro.c bezctx.c spiroentrypoints.c repro/regress.c \
   -lm -o repro/regress_new
./repro/regress_new > /tmp/new.txt && diff /tmp/orig.txt /tmp/new.txt
```

官方测试：`autoreconf -i && ./configure && make check`（24/24 通过，
`--enable-test_inputs` 构建同样 24/24）。

## 4. 数据文件

- `scan-results-original.txt` / `scan-results-fixed.txt`：456 用例边界扫描前后。
- `fuzz-results-original.txt` / `fuzz-results-fixed.txt`：10 万次随机 fuzz 前后。

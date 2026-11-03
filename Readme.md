# SABRE: Influence Minimization via Edge Blocking

## Build

Requires a C++17 compiler.

```bash
make        # produces ./IMIN_EB
```

## Usage

```bash
./IMIN_EB -dataset <path>/ -rumorNum <id> -k <#edges> -algo <name> [options]
```

Examples (on the included toy dataset):

```bash
./IMIN_EB -dataset datasets/toy/ -rumorNum 5 -k 5 -algo SABRE -epsilon 0.1 -gamma 0.3
./IMIN_EB -dataset datasets/toy/ -rumorNum 5 -k 5 -algo Gre-DT+
./IMIN_EB -dataset datasets/toy/ -rumorNum 5 -k 5 -algo "Gre*-DT"
```

| Argument      | Description                                                                           |
| ------------- | ------------------------------------------------------------------------------------- |
| `-dataset`  | Dataset folder path                                                                   |
| `-rumorNum` | Seed set file ID (loads `rumorSet_<id>.txt`)                                        |
| `-k`        | Number of edges to block                                                              |
| `-algo`     | `SABRE`, `Gre-MC`, `Gre-DT`, `Gre*-DT`, `Gre-DT+`, `ProbDeg`, `SeedAdj` |
| `-epsilon`  | SABRE: error tolerance ε (default 0.1)                                               |
| `-gamma`    | SABRE: objective-estimator relative error γ (default 0.3)                            |
| `-delta`    | SABRE: failure probability δ (default 1/n)                                           |
| `-theta`    | Gre-MC / Gre-DT / Gre*-DT / Gre-DT+: sample size (default 10000)                      |
| `-res`      | Result output file path (optional)                                                    |

## Dataset Format

Each dataset folder contains (see `datasets/toy/` for an example):

- `attribute.txt` — two lines: `n=<#nodes>` and `m=<#edges>`
- `graph.txt` — one directed edge per line: `source target probability` (an optional `n m` header line is auto-detected and skipped)
- `rumorSet_<id>.txt` — seed node IDs, one per line

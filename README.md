# LESS: Linear Equations Systems Solver

Academic activity for the parallel programming lab course from the Computer Science undegraduate course of Universidade Federal Fluminense.

## How tos

To check for syntax errors, you may use:

`make check`

To generate a release build (placed at `build/release/main`), just use:

`make release`

To generate a release build and run, use:

`make run p={PROCESS_NUMBERS}`

## Results

| Matrix size | Time - 1 process | Time - 2 processes | Speedup - 2 processes | Time - 4 processes | Speedup - 4 processes |
|:-----------:|------------------|--------------------|-----------------------|--------------------|-----------------------|
| 50          | 0,294 s          | 0,280 s            | 1,05                  | 0,306              | 0,96                  |
| 500         | 0,542 s          | 0,450 s            | 1,204                 | 0,394 s            | 1,375                 |
| 1500        | 5,814 s          | 3,188 s            | 1,823                 | 1,751 s            | 3,320                 |
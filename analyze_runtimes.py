import re
import glob
from collections import defaultdict
from datetime import timedelta
import statistics

runtimes = defaultdict(list)

for filename in glob.glob("output_*.txt"):
    with open(filename, "r", errors="ignore") as f:
        content = f.read()

    # Skip unfinished jobs
    runtime_match = re.search(
        r"Total executation took\s+(\d{2}):(\d{2}):(\d{2})",
        content
    )

    if not runtime_match:
        continue

    # Extract L
    L_match = re.search(r"/L=([\d.]+)/", content)

    if not L_match:
        continue

    L = float(L_match.group(1))

    # Convert runtime to seconds
    hours, minutes, seconds = map(int, runtime_match.groups())
    runtime = hours * 3600 + minutes * 60 + seconds

    runtimes[L].append(runtime)


# Print results
print(f"{'L':>8} {'N':>6} {'Mean [s]':>14} {'Std [s]':>14} {'Mean [h:m:s]':>16}")

for L in sorted(runtimes):
    values = runtimes[L]
    mean = statistics.mean(values)
    std = statistics.stdev(values) if len(values) > 1 else 0.0

    mean_formatted = str(timedelta(seconds=round(mean)))

    print(
        f"{L:8g} {len(values):6d} "
        f"{mean:14.1f} {std:14.1f} "
        f"{mean_formatted:>16}"
    )

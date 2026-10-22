#!/usr/bin/env python3
import csv
import statistics
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

base = Path(__file__).resolve().parent
src = base / 'results' / 'raw.csv'
rows = list(csv.DictReader(src.open(encoding='utf8', newline='')))
groups = defaultdict(list)
for row in rows:
    key = (row['experiment'], row['variant'], int(row['size']), int(row['threads']), row['unit'])
    groups[key].append(float(row['metric']))
summary = base / 'results' / 'summary_median.csv'
with summary.open('w', newline='', encoding='utf8') as file:
    writer = csv.writer(file)
    writer.writerow(['experiment','variant','size','threads','median_metric','unit','samples'])
    for key, values in sorted(groups.items()):
        writer.writerow([*key[:4], f'{statistics.median(values):.6f}', key[4],len(values)])

def plot(exp, x_label, y_label, output, log_x=False):
    data=defaultdict(list)
    for (experiment,variant,size,threads,unit), vals in groups.items():
        if experiment!=exp: continue
        x=threads if exp in ('false_sharing','sync') else size
        if exp=='stride': x=int(variant.split('_')[1])
        data[variant if exp not in ('working_set','matrix') else variant].append((x,statistics.median(vals)))
    fig, ax=plt.subplots(figsize=(9,5.5))
    for name, points in sorted(data.items()):
        points.sort()
        ax.plot([p[0] for p in points],[p[1] for p in points], marker='o',label=name)
    ax.set_xlabel(x_label);ax.set_ylabel(y_label);ax.set_title(exp.replace('_',' ').title())
    if log_x:ax.set_xscale('log',base=2)
    ax.grid(alpha=.3);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(base/'plots'/output,dpi=160);plt.close(fig)

(base/'plots').mkdir(exist_ok=True)
plot('working_set','Working set (bytes)','ns per access','01_working_set.png',True)
plot('stride','Stride (bytes)','ns per access','02_stride.png',True)
plot('matrix','Matrix dimension N','GFLOP/s','03_matrix.png')
plot('false_sharing','Threads','Million atomic increments/s','04_false_sharing.png')
plot('sync','Threads','Correct updates (%)','05_sync.png')
plot('volatile','Iterations','ns per increment','06_volatile.png')
print(f'Saved {summary} and six figures in plots/')

#!/usr/bin/env python3
"""
Generate comprehensive HTML report for genome assembly quality assessment
"""

import argparse
import base64
import io
import re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def parse_assembly_stats(stats_file):
    """Parse assembly stats file"""
    df = pd.read_csv(stats_file, sep='\t')
    return df

def parse_busco_summary(busco_file):
    """Parse BUSCO short_summary file"""
    with open(busco_file, 'r') as f:
        content = f.read()

    # Extract sample name from filename
    sample_name = Path(busco_file).stem.replace('short_summary.specific.mammalia_odb12.', '').replace('.busco', '')

    # Parse the summary line
    summary_pattern = r'C:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%'
    match = re.search(summary_pattern, content)

    if match:
        complete_pct = float(match.group(1))
        single_pct = float(match.group(2))
        duplicated_pct = float(match.group(3))
        fragmented_pct = float(match.group(4))
        missing_pct = float(match.group(5))
    else:
        return None

    # Extract counts
    complete_match = re.search(r'(\d+)\s+Complete BUSCOs', content)
    single_match = re.search(r'(\d+)\s+Complete and single-copy BUSCOs', content)
    duplicated_match = re.search(r'(\d+)\s+Complete and duplicated BUSCOs', content)
    fragmented_match = re.search(r'(\d+)\s+Fragmented BUSCOs', content)
    missing_match = re.search(r'(\d+)\s+Missing BUSCOs', content)
    total_match = re.search(r'(\d+)\s+Total BUSCO groups searched', content)

    return {
        'sample': sample_name,
        'complete_pct': complete_pct,
        'single_pct': single_pct,
        'duplicated_pct': duplicated_pct,
        'fragmented_pct': fragmented_pct,
        'missing_pct': missing_pct,
        'complete': int(complete_match.group(1)) if complete_match else 0,
        'single': int(single_match.group(1)) if single_match else 0,
        'duplicated': int(duplicated_match.group(1)) if duplicated_match else 0,
        'fragmented': int(fragmented_match.group(1)) if fragmented_match else 0,
        'missing': int(missing_match.group(1)) if missing_match else 0,
        'total': int(total_match.group(1)) if total_match else 0
    }

def parse_repeatmasker_tbl(tbl_file):
    """Parse RepeatMasker .tbl file"""
    with open(tbl_file, 'r') as f:
        content = f.read()

    # Extract sample name from filename
    sample_name = Path(tbl_file).stem.replace('.repeats', '')

    # Extract total length
    total_length_match = re.search(r'total length:\s+(\d+) bp', content)
    total_length = int(total_length_match.group(1)) if total_length_match else 0

    # Extract bases masked
    masked_match = re.search(r'bases masked:\s+(\d+) bp', content)
    bases_masked = int(masked_match.group(1)) if masked_match else 0

    # Extract repeat categories
    repeats = {}

    # Retroelements
    retro_match = re.search(r'Retroelements\s+\d+\s+(\d+) bp', content)
    if retro_match:
        repeats['Retroelements'] = int(retro_match.group(1))

    # DNA transposons
    dna_match = re.search(r'DNA transposons\s+\d+\s+(\d+) bp', content)
    if dna_match:
        repeats['DNA transposons'] = int(dna_match.group(1))

    # Rolling-circles
    rolling_match = re.search(r'Rolling-circles\s+\d+\s+(\d+) bp', content)
    if rolling_match:
        repeats['Rolling-circles'] = int(rolling_match.group(1))

    # Unclassified
    unclass_match = re.search(r'Unclassified:\s+\d+\s+(\d+) bp', content)
    if unclass_match:
        repeats['Unclassified'] = int(unclass_match.group(1))

    # Simple repeats
    simple_match = re.search(r'Simple repeats:\s+\d+\s+(\d+) bp', content)
    if simple_match:
        repeats['Simple repeats'] = int(simple_match.group(1))

    # Low complexity
    low_match = re.search(r'Low complexity:\s+\d+\s+(\d+) bp', content)
    if low_match:
        repeats['Low complexity'] = int(low_match.group(1))

    # Calculate non-repetitive
    total_repeats = sum(repeats.values())
    repeats['Non-repetitive'] = total_length - total_repeats

    return {
        'sample': sample_name,
        'total_length': total_length,
        'bases_masked': bases_masked,
        'repeats': repeats
    }

def parse_merqury_qv(qv_file):
    """Parse merqury QV file to extract genome-wide QV value"""
    with open(qv_file, 'r') as f:
        lines = f.readlines()

    # QV file format: assembly_name num_kmers error_rate qv
    # We want the overall QV (first line usually has the total)
    if lines:
        parts = lines[0].strip().split()
        if len(parts) >= 4:
            return {
                'assembly': parts[0],
                'num_kmers': int(parts[1]),
                'error_rate': float(parts[2]),
                'qv': float(parts[3])
            }
    return None

def plot_assembly_stats(stats_files, sample_name):
    """Create assembly statistics plots for a single sample with multiple assemblies"""
    all_stats = []
    for f in stats_files:
        df = parse_assembly_stats(f)
        # Extract assembly type from filename (hap1, hap2, or p_ctg)
        filename = Path(f).stem.replace('.stats', '')
        if 'hap1' in filename:
            assembly_type = 'Hap1'
        elif 'hap2' in filename:
            assembly_type = 'Hap2'
        else:
            assembly_type = 'Primary'
        df['assembly'] = assembly_type
        all_stats.append(df)

    combined = pd.concat(all_stats, ignore_index=True)

    # Create subplots for each metric
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'Assembly Statistics - {sample_name}', fontsize=16, fontweight='bold')

    metrics = ['total_length', 'longest_sequence', 'num_sequence', 'N50', 'L50']

    # Get unique assemblies for x-axis in specific order: Primary, Hap1, Hap2
    unique_assemblies = combined['assembly'].unique()
    assembly_order = ['Primary', 'Hap1', 'Hap2']
    assemblies = [asm for asm in assembly_order if asm in unique_assemblies]
    x = np.arange(len(assemblies))
    width = 0.35

    for idx, metric in enumerate(metrics):
        row = idx // 3
        col = idx % 3
        ax = axes[row, col]

        # Separate contigs and scaffolds
        contigs_data = []
        scaffolds_data = []
        for asm in assemblies:
            asm_data = combined[combined['assembly'] == asm]
            contig_val = asm_data[asm_data['level'] == 'contigs'][metric].values
            scaffold_val = asm_data[asm_data['level'] == 'scaffolds'][metric].values
            contigs_data.append(contig_val[0] if len(contig_val) > 0 else 0)
            scaffolds_data.append(scaffold_val[0] if len(scaffold_val) > 0 else 0)

        ax.bar(x - width/2, contigs_data, width, label='Contigs', alpha=0.8)
        ax.bar(x + width/2, scaffolds_data, width, label='Scaffolds', alpha=0.8)

        ax.set_xlabel('Assembly')
        ax.set_ylabel(metric.replace('_', ' ').title())
        ax.set_title(metric.replace('_', ' ').title())
        ax.set_xticks(x)
        ax.set_xticklabels(assemblies, rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3)

    # Remove the extra subplot
    fig.delaxes(axes[1, 2])

    plt.tight_layout()

    # Convert to base64
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close()

    return img_base64

def plot_busco_results(busco_files, sample_name):
    """Create BUSCO plot for a single sample with multiple assemblies"""
    busco_data = []
    for f in busco_files:
        data = parse_busco_summary(f)
        if data:
            # Extract assembly type from sample name
            if 'hap1' in data['sample']:
                data['assembly'] = 'Hap1'
            elif 'hap2' in data['sample']:
                data['assembly'] = 'Hap2'
            else:
                data['assembly'] = 'Primary'
            busco_data.append(data)

    if not busco_data:
        return None

    df = pd.DataFrame(busco_data)

    # Order assemblies: Primary, Hap1, Hap2
    assembly_order = ['Primary', 'Hap1', 'Hap2']
    df['assembly'] = pd.Categorical(df['assembly'], categories=assembly_order, ordered=True)
    df = df.sort_values('assembly')

    fig, ax = plt.subplots(figsize=(max(10, len(df) * 1.5), 8))

    x = np.arange(len(df))
    width = 0.6

    # Plot stacked bars
    p1 = ax.bar(x, df['single_pct'], width, label='Complete (Single)', color='#56B4E9')
    p2 = ax.bar(x, df['duplicated_pct'], width, bottom=df['single_pct'],
                label='Complete (Duplicated)', color='#3274A1')
    p3 = ax.bar(x, df['fragmented_pct'], width,
                bottom=df['single_pct'] + df['duplicated_pct'],
                label='Fragmented', color='#F0E442')
    p4 = ax.bar(x, df['missing_pct'], width,
                bottom=df['single_pct'] + df['duplicated_pct'] + df['fragmented_pct'],
                label='Missing', color='#E69F00')

    ax.set_ylabel('Percentage (%)', fontsize=12)
    ax.set_title(f'BUSCO Completeness Assessment - {sample_name}', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(df['assembly'], rotation=45, ha='right')
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_ylim(0, 100)
    ax.grid(True, axis='y', alpha=0.3)

    # Add percentage labels (complete single-copy percentage)
    for pos, (_, row) in enumerate(df.iterrows()):
        if row['single_pct'] > 5:
            ax.text(pos, row['single_pct']/2, f"{row['single_pct']:.1f}%",
                   ha='center', va='center', fontweight='bold', fontsize=10)

    plt.tight_layout()

    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close()

    return img_base64

def plot_repeat_content(repeat_files, sample_name):
    """Create repeat content pie charts for a single sample with multiple assemblies"""
    repeat_data = []
    for f in repeat_files:
        data = parse_repeatmasker_tbl(f)
        # Extract assembly type from filename
        if 'hap1' in data['sample']:
            data['assembly'] = 'Hap1'
        elif 'hap2' in data['sample']:
            data['assembly'] = 'Hap2'
        else:
            data['assembly'] = 'Primary'
        repeat_data.append(data)

    if not repeat_data:
        return None

    n_assemblies = len(repeat_data)
    cols = min(3, n_assemblies)
    rows = (n_assemblies + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 6*rows))
    if n_assemblies == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if n_assemblies > 1 else [axes]

    colors = plt.cm.Set3(np.linspace(0, 1, 10))

    for idx, data in enumerate(repeat_data):
        ax = axes[idx]

        repeats = data['repeats']
        # Only include categories with > 0 bp
        filtered_repeats = {k: v for k, v in repeats.items() if v > 0}

        labels = list(filtered_repeats.keys())
        sizes = list(filtered_repeats.values())

        # Create pie chart without labels (will use legend instead)
        wedges, texts, autotexts = ax.pie(sizes, autopct='%1.1f%%',
                                            colors=colors[:len(labels)], startangle=90)

        # Make percentage text more readable
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
            autotext.set_fontsize(9)

        # Add colored legend
        ax.legend(wedges, labels, loc='center left', bbox_to_anchor=(1, 0, 0.5, 1),
                 frameon=True, fontsize=9)

        ax.set_title(f"{data['assembly']}\nTotal: {data['total_length']:,} bp",
                    fontweight='bold')

    # Remove extra subplots
    for idx in range(n_assemblies, len(axes)):
        fig.delaxes(axes[idx])

    plt.suptitle(f'Repeat Content Analysis - {sample_name}', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()

    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close()

    return img_base64

def encode_image(image_path):
    """Encode image file to base64"""
    with open(image_path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')

def generate_html_report(output_file, sample_name, stats_files, busco_files, merqury_files, merqury_qv_files, repeat_files):
    """Generate comprehensive HTML report for a single sample"""

    # Generate plots
    stats_plot = plot_assembly_stats(stats_files, sample_name)
    busco_plot = plot_busco_results(busco_files, sample_name)
    repeat_plot = plot_repeat_content(repeat_files, sample_name)

    # Parse QV values
    qv_data = {}
    for qv_file in merqury_qv_files:
        qv_info = parse_merqury_qv(qv_file)
        if qv_info:
            filename = Path(qv_file).stem
            if 'diploid' in filename:
                qv_type = 'Diploid'
            else:
                qv_type = 'Primary'
            qv_data[qv_type] = qv_info

    # Encode merqury images - only diploid spectra-asm plot
    diploid_asm_image = None
    for mf in merqury_files:
        if mf.endswith('.png') and 'diploid' in str(mf) and 'spectra-asm' in str(mf):
            diploid_asm_image = encode_image(mf)
            break

    # Create HTML
    html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Genome Assembly Quality Report - {sample_name}</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            box-shadow: 0 0 10px rgba(0,0,0,0.1);
            border-radius: 8px;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 40px;
            border-bottom: 2px solid #bdc3c7;
            padding-bottom: 8px;
        }}
        .plot-container {{
            margin: 30px 0;
            text-align: center;
        }}
        .plot-container img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
            padding: 5px;
        }}
        .merqury-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .merqury-item {{
            text-align: center;
        }}
        .merqury-item h3 {{
            color: #34495e;
            margin-bottom: 10px;
        }}
        .qv-info {{
            background-color: #e8f4f8;
            border-left: 4px solid #3498db;
            padding: 15px;
            margin: 20px 0;
            border-radius: 4px;
        }}
        .qv-info strong {{
            color: #2c3e50;
        }}
        .timestamp {{
            color: #7f8c8d;
            font-size: 0.9em;
            margin-top: 40px;
            text-align: center;
            padding-top: 20px;
            border-top: 1px solid #ecf0f1;
        }}
        .section-description {{
            color: #555;
            font-style: italic;
            margin-bottom: 20px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Genome Assembly Quality Report - {sample_name}</h1>

        <h2>📊 Assembly Statistics</h2>
        <p class="section-description">
            Comparison of key assembly metrics across the three assemblies, showing both contig and scaffold-level statistics.
        </p>
        <div class="plot-container">
            <img src="data:image/png;base64,{stats_plot}" alt="Assembly Statistics">
        </div>

        <h2>🎯 BUSCO Completeness</h2>
        <p class="section-description">
            Benchmarking Universal Single-Copy Orthologs (BUSCO) assessment shows the completeness and quality of each assembly.
        </p>
        <div class="plot-container">
            <img src="data:image/png;base64,{busco_plot}" alt="BUSCO Results">
        </div>

        <h2>📈 Merqury K-mer Analysis</h2>
        <p class="section-description">
            K-mer spectra plots show the copy number distribution in the assembly, helping assess assembly accuracy and completeness.
        </p>
"""

    # Add QV values section
    if qv_data:
        html += """        <div class="qv-info">
            <strong>Genome-wide Quality Values (QV):</strong><br>
"""
        for qv_type, qv_info in sorted(qv_data.items()):
            html += f"""            • {qv_type}: QV = {qv_info['qv']:.2f} (Error rate: {qv_info['error_rate']:.2e})<br>
"""
        html += """        </div>
"""

    # Add diploid merqury assembly spectra plot
    if diploid_asm_image:
        html += f"""
        <div class="plot-container">
            <img src="data:image/png;base64,{diploid_asm_image}" alt="Diploid assembly k-mer spectra" style="max-width: 800px;">
        </div>
"""
    else:
        html += """
        <p style="color: #e74c3c; font-style: italic;">Diploid assembly spectra plot not available.</p>
"""

    html += f"""

        <h2>🔄 Repeat Content</h2>
        <p class="section-description">
            RepeatMasker analysis showing the distribution of repetitive elements and non-repetitive regions in the genome.
        </p>
        <div class="plot-container">
            <img src="data:image/png;base64,{repeat_plot}" alt="Repeat Content">
        </div>

        <div class="timestamp">
            Report generated on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
        </div>
    </div>
</body>
</html>
"""

    with open(output_file, 'w') as f:
        f.write(html)

    print(f"Report generated: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate assembly quality report for a single sample')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--stats', nargs='+', required=True, help='Assembly stats files')
    parser.add_argument('--busco', nargs='+', required=True, help='BUSCO summary files')
    parser.add_argument('--merqury-png', nargs='+', required=True, help='Merqury spectra PNG files')
    parser.add_argument('--merqury-qv', nargs='+', required=True, help='Merqury QV files')
    parser.add_argument('--repeats', nargs='+', required=True, help='RepeatMasker .tbl files')
    parser.add_argument('--output', required=True, help='Output HTML file')

    args = parser.parse_args()

    generate_html_report(args.output, args.sample, args.stats, args.busco,
                        args.merqury_png, args.merqury_qv, args.repeats)

if __name__ == '__main__':
    main()

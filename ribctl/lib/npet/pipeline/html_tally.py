import os
import json
import argparse
from pathlib import Path
from collections import Counter

def generate_processing_report(log_dir):
    """Generate HTML report for NPET mesh processing logs"""
    log_dir = Path(log_dir)

    # Find all processing logs
    log_files = list(log_dir.glob("*_processing_log.json"))

    if not log_files:
        print(f"No processing logs found in {log_dir}")
        return

    print(f"Found {len(log_files)} processing logs")

    # Load all logs
    structures = []
    for log_file in log_files:
        try:
            with open(log_file, 'r') as f:
                data = json.load(f)
                structures.append(data)
        except Exception as e:
            print(f"Error loading {log_file}: {e}")

    # Count successes and failures
    succeeded = [s for s in structures if s["overall_status"] == "success"]
    failed = [s for s in structures if s["overall_status"] != "success"]

    # Count failures by stage
    failed_stages = Counter()
    for struct in failed:
        failed_stage = struct["summary"].get("failed_stage")
        if failed_stage:
            failed_stages[failed_stage] += 1

    # Generate HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>NPET Mesh Processing Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }}
        h1, h2, h3 {{
            color: #2c3e50;
        }}
        .summary {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .card {{
            background-color: #fff;
            border-radius: 8px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
            padding: 20px;
            text-align: center;
        }}
        .success {{
            background-color: #d4edda;
            color: #155724;
        }}
        .failed {{
            background-color: #f8d7da;
            color: #721c24;
        }}
        .card h3 {{
            margin-top: 0;
            font-size: 18px;
        }}
        .card .number {{
            font-size: 48px;
            font-weight: bold;
            margin: 10px 0;
            line-height: 1;
        }}
        .card .percent {{
            font-size: 18px;
            opacity: 0.8;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 30px;
        }}
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #f8f9fa;
            font-weight: 600;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .status-badge {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 14px;
        }}
        .status-success {{
            background-color: #d4edda;
            color: #155724;
        }}
        .status-failure {{
            background-color: #f8d7da;
            color: #721c24;
        }}
        .charts {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            margin-bottom: 30px;
        }}
        .chart-container {{
            height: 400px;
        }}
        .error-message {{
            max-width: 500px;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }}
        .error-message:hover {{
            white-space: normal;
            overflow: visible;
        }}
        .toggle-section {{
            margin-bottom: 20px;
        }}
        .toggle-button {{
            background-color: #4CAF50;
            border: none;
            color: white;
            padding: 10px 15px;
            text-align: center;
            text-decoration: none;
            display: inline-block;
            font-size: 16px;
            margin: 4px 2px;
            cursor: pointer;
            border-radius: 4px;
        }}
        #failuresTable {{
            display: none;
        }}
        #successTable {{
            display: none;
        }}
    </style>
</head>
<body>
    <h1>NPET Mesh Processing Report</h1>

    <div class="summary">
        <div class="card">
            <h3>Total Structures</h3>
            <div class="number">{len(structures)}</div>
        </div>
        <div class="card success">
            <h3>Success</h3>
            <div class="number">{len(succeeded)}</div>
            <div class="percent">{len(succeeded)/len(structures)*100:.1f}%</div>
        </div>
        <div class="card failed">
            <h3>Failed</h3>
            <div class="number">{len(failed)}</div>
            <div class="percent">{len(failed)/len(structures)*100:.1f}%</div>
        </div>
    </div>

    <div class="charts">
        <div class="chart-container">
            <canvas id="stageFailureChart"></canvas>
        </div>
        <div class="chart-container">
            <canvas id="completionChart"></canvas>
        </div>
    </div>

    <div class="toggle-section">
        <button class="toggle-button" onclick="toggleTable('failuresTable')">Show/Hide Failures</button>
        <button class="toggle-button" onclick="toggleTable('successTable')">Show/Hide Successes</button>
    </div>

    <h2>Failed Structures ({len(failed)})</h2>
    <table id="failuresTable">
        <thead>
            <tr>
                <th>RCSB ID</th>
                <th>Failed Stage</th>
                <th>Duration (s)</th>
                <th>Error</th>
            </tr>
        </thead>
        <tbody>
"""

    # Add failed structures to table
    for struct in sorted(failed, key=lambda x: x["rcsb_id"]):
        rcsb_id = struct["rcsb_id"]
        failed_stage = struct["summary"].get("failed_stage", "Unknown")
        duration = struct.get("duration", 0)

        # Find the error message
        error_msg = "Unknown error"
        if failed_stage and failed_stage in struct["stages"]:
            stage_info = struct["stages"][failed_stage]
            if stage_info.get("error") and stage_info["error"].get("message"):
                error_msg = stage_info["error"]["message"]

        html += f"""
            <tr>
                <td>{rcsb_id}</td>
                <td>{failed_stage}</td>
                <td>{duration}</td>
                <td class="error-message">{error_msg}</td>
            </tr>"""

    html += """
        </tbody>
    </table>

    <h2>Successful Structures ({len(succeeded)})</h2>
    <table id="successTable">
        <thead>
            <tr>
                <th>RCSB ID</th>
                <th>Duration (s)</th>
                <th>Watertight</th>
                <th>Artifacts</th>
            </tr>
        </thead>
        <tbody>
"""

    # Add successful structures to table
    for struct in sorted(succeeded, key=lambda x: x["rcsb_id"]):
        rcsb_id = struct["rcsb_id"]
        duration = struct.get("duration", 0)
        watertight = "Yes" if struct["summary"].get("watertight", False) else "No"
        artifact_count = len(struct["summary"].get("artifacts_generated", []))

        html += f"""
            <tr>
                <td>{rcsb_id}</td>
                <td>{duration:.2f}</td>
                <td>{watertight}</td>
                <td>{artifact_count}</td>
            </tr>"""

    # Prepare data for charts
    stages_data = []
    for stage, count in failed_stages.items():
        stages_data.append({"stage": stage, "count": count})

    # Complete the HTML
    html += """
        </tbody>
    </table>

    <script>
        // Failed stages chart
        const stageCtx = document.getElementById('stageFailureChart').getContext('2d');
        const stageFailureChart = new Chart(stageCtx, {
            type: 'bar',
            data: {
                labels: [""" + ', '.join([f'"{stage["stage"]}"' for stage in stages_data]) + """],
                datasets: [{
                    label: 'Failed Structures by Stage',
                    data: [""" + ', '.join([str(stage["count"]) for stage in stages_data]) + """],
                    backgroundColor: 'rgba(255, 99, 132, 0.5)',
                    borderColor: 'rgba(255, 99, 132, 1)',
                    borderWidth: 1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Number of Failures'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Processing Stage'
                        }
                    }
                }
            }
        });

        // Completion status chart
        const completionCtx = document.getElementById('completionChart').getContext('2d');
        const completionChart = new Chart(completionCtx, {
            type: 'pie',
            data: {
                labels: ['Success', 'Failed'],
                datasets: [{
                    data: [""" + f"{len(succeeded)}, {len(failed)}" + """],
                    backgroundColor: [
                        'rgba(75, 192, 192, 0.5)',
                        'rgba(255, 99, 132, 0.5)'
                    ],
                    borderColor: [
                        'rgba(75, 192, 192, 1)',
                        'rgba(255, 99, 132, 1)'
                    ],
                    borderWidth: 1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: {
                        position: 'top',
                    },
                    title: {
                        display: true,
                        text: 'Structure Processing Status'
                    }
                }
            }
        });

        function toggleTable(tableId) {
            const table = document.getElementById(tableId);
            if (table.style.display === "none" || table.style.display === "") {
                table.style.display = "table";
            } else {
                table.style.display = "none";
            }
        }
    </script>
</body>
</html>
"""

    # Write HTML to file
    output_path = log_dir / "npet_processing_report.html"
    with open(output_path, 'w') as f:
        f.write(html)

    print(f"Report generated at {output_path}")
    return output_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate NPET mesh processing report")
    parser.add_argument("--log-dir", default="logs", help="Directory containing processing logs")
    args = parser.parse_args()

    generate_processing_report(args.log_dir)

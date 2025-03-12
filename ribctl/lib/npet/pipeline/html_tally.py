#!/usr/bin/env python3
"""
Simplified Ribosome Pipeline Dashboard Generator

This script generates an HTML dashboard from ribosome pipeline processing logs.
With specific focus on column sorting and structure selection.
"""

import os
import sys
import json
import datetime
from pathlib import Path
from typing import Dict, List, Any

# Default stages in expected processing order
DEFAULT_STAGES = [
    "setup", 
    "ptc_identification", 
    "constriction_identification",
    "alpha_shape", 
    "landmark_identification",
    "entity_filtering", 
    "point_cloud_processing",
    "clustering", 
    "refinement", 
    "surface_extraction",
    "normal_estimation", 
    "mesh_reconstruction",
    "validation", 
    "complete"
]

def read_log_file(file_path: Path) -> Dict[str, Any]:
    """Read and parse a JSON log file."""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return {}

def scan_logs_directory(directory: Path) -> List[Dict[str, Any]]:
    """Scan directory for processing logs and return parsed data."""
    log_data = []
    
    # Find all JSON files in the directory
    log_files = list(directory.glob("*_processing_log.json"))
    
    if not log_files:
        # If no *_processing_log.json files, try all JSON files
        log_files = list(directory.glob("*.json"))
    
    for file_path in log_files:
        data = read_log_file(file_path)
        if data and 'rcsb_id' in data and 'stages' in data:
            log_data.append(data)
    
    # Sort by RCSB ID
    log_data.sort(key=lambda x: x.get('rcsb_id', ''))
    return log_data

def format_duration(seconds: float) -> str:
    """Format duration in seconds to human-readable string."""
    if seconds is None:
        return "N/A"
    
    if seconds < 0.001:
        return "<0.001s"
    elif seconds < 1:
        return f"{seconds:.3f}s"
    elif seconds < 60:
        return f"{seconds:.2f}s"
    elif seconds < 3600:
        minutes = int(seconds // 60)
        remaining_seconds = seconds % 60
        return f"{minutes}m {remaining_seconds:.1f}s"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        return f"{hours}h {minutes}m"

def generate_html(logs: List[Dict[str, Any]], output_path: Path):
    """Generate HTML dashboard from log data."""
    if not logs:
        print("No valid log data found.")
        return
    
    # Get all unique stage names from logs
    all_stages = set()
    for log in logs:
        all_stages.update(log.get('stages', {}).keys())
    
    # Order stages according to DEFAULT_STAGES with any additional stages at the end
    stages = [stage for stage in DEFAULT_STAGES if stage in all_stages]
    stages.extend(sorted(all_stages - set(stages)))
    
    # Count success and failure
    success_count = sum(1 for log in logs if log.get('overall_status') == 'success')
    failure_count = sum(1 for log in logs if log.get('overall_status') == 'failure')
    
    # Calculate average duration for successful runs
    successful_durations = [log.get('duration', 0) for log in logs if log.get('overall_status') == 'success']
    avg_duration = sum(successful_durations) / len(successful_durations) if successful_durations else 0
    total_duration = sum(successful_durations) if successful_durations else 0
    
    # Find fastest and slowest
    if successful_durations:
        fastest_id = min((log for log in logs if log.get('overall_status') == 'success'), 
                       key=lambda x: x.get('duration', float('inf')))['rcsb_id']
        fastest_time = min(successful_durations)
        
        slowest_id = max((log for log in logs if log.get('overall_status') == 'success'), 
                       key=lambda x: x.get('duration', 0))['rcsb_id']
        slowest_time = max(successful_durations)
    else:
        fastest_id, fastest_time = "N/A", 0
        slowest_id, slowest_time = "N/A", 0
    
    # Generate HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Ribosome Pipeline Processing Dashboard</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans', 'Helvetica Neue', sans-serif;
            line-height: 1.6;
            color: #333;
            margin: 0;
            padding: 20px;
            background-color: #f5f7fa;
        }}
        
        h1 {{
            color: #2c3e50;
            border-bottom: 2px solid #3498db;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        
        .controls {{
            margin-bottom: 20px;
            display: flex;
            gap: 15px;
            align-items: center;
            flex-wrap: wrap;
        }}
        
        .search-box {{
            padding: 8px 12px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
            min-width: 200px;
        }}
        
        .filter-select {{
            padding: 8px 12px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }}
        
        .dashboard {{
            overflow-x: auto;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.05);
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 14px;
        }}
        
        th {{
            background-color: #f8f9fa;
            text-align: left;
            padding: 12px;
            position: sticky;
            top: 0;
            z-index: 1;
            border-bottom: 2px solid #ddd;
            cursor: pointer;
        }}
        
        th.sorted-asc:after {{
            content: " ↑";
            color: #3498db;
        }}
        
        th.sorted-desc:after {{
            content: " ↓";
            color: #3498db;
        }}
        
        td {{
            padding: 10px 12px;
            border-bottom: 1px solid #eee;
            white-space: nowrap;
        }}
        
        tr:hover {{
            background-color: #f5f5f5;
        }}
        
        tr.selected {{
            background-color: #e3f2fd !important;
        }}
        
        .status-cell {{
            width: 30px;
            text-align: center;
        }}
        
        .success {{
            background-color: #d4edda;
            color: #155724;
        }}
        
        .failure {{
            background-color: #f8d7da;
            color: #721c24;
        }}
        
        .pending {{
            background-color: #fff3cd;
            color: #856404;
        }}
        
        .skipped {{
            background-color: #e2e3e5;
            color: #383d41;
        }}
        
        .in-progress {{
            background-color: #cce5ff;
            color: #004085;
        }}
        
        .status-icon {{
            font-size: 14px;
            font-weight: bold;
        }}
        
        .duration {{
            color: #6c757d;
            font-size: 12px;
        }}
        
        .summary {{
            display: flex;
            gap: 20px;
            margin-bottom: 20px;
            flex-wrap: wrap;
        }}
        
        .summary-card {{
            background: white;
            border-radius: 8px;
            padding: 15px;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.05);
            flex: 1;
            min-width: 200px;
        }}
        
        .summary-title {{
            font-size: 14px;
            color: #6c757d;
            margin-bottom: 8px;
        }}
        
        .summary-value {{
            font-size: 28px;
            font-weight: bold;
            color: #2c3e50;
        }}
        
        .summary-subvalue {{
            font-size: 14px;
            color: #6c757d;
            margin-top: 5px;
        }}
        
        .footer {{
            text-align: center;
            margin-top: 40px;
            color: #6c757d;
            font-size: 14px;
        }}
        
        /* Selection panel styles */
        #selection-panel {{
            position: fixed;
            bottom: 0;
            left: 0;
            right: 0;
            background: #fff;
            padding: 15px;
            border-top: 1px solid #ddd;
            box-shadow: 0 -2px 10px rgba(0,0,0,0.1);
            display: none;
            z-index: 1000;
        }}
        
        #selection-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 10px;
        }}
        
        #selection-count {{
            background: #3498db;
            color: white;
            border-radius: 50%;
            width: 24px;
            height: 24px;
            display: inline-flex;
            align-items: center;
            justify-content: center;
            margin-left: 10px;
        }}
        
        #selection-text {{
            border: 1px solid #ddd;
            padding: 10px;
            max-height: 200px;
            overflow-y: auto;
            font-family: monospace;
            background: #f8f9fa;
            white-space: pre;
        }}
        
        .btn {{
            padding: 8px 15px;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            margin-left: 10px;
        }}
        
        .btn-primary {{
            background: #3498db;
            color: white;
        }}
        
        .btn-secondary {{
            background: #e2e3e5;
            color: #333;
        }}
    </style>
</head>
<body>
    <h1>Ribosome Pipeline Processing Dashboard</h1>
    
    <div class="summary">
        <div class="summary-card">
            <div class="summary-title">Total Structures</div>
            <div class="summary-value">{len(logs)}</div>
        </div>
        <div class="summary-card">
            <div class="summary-title">Success Rate</div>
            <div class="summary-value">{success_count / len(logs) * 100:.1f}%</div>
            <div class="summary-subvalue">{success_count} successful, {failure_count} failed</div>
        </div>
        <div class="summary-card">
            <div class="summary-title">Average Processing Time</div>
            <div class="summary-value">{format_duration(avg_duration)}</div>
            <div class="summary-subvalue">Total: {format_duration(total_duration)}</div>
        </div>
        <div class="summary-card">
            <div class="summary-title">Processing Stats</div>
            <div class="summary-subvalue">
                <div>Fastest: {fastest_id} ({format_duration(fastest_time)})</div>
                <div>Slowest: {slowest_id} ({format_duration(slowest_time)})</div>
            </div>
        </div>
    </div>

    <div class="controls">
        <input type="text" id="search-input" class="search-box" placeholder="Search by ID...">
        <button onclick="clearSelectedRows()" class="btn btn-secondary">Clear Selection</button>
    </div>
    
    <div class="dashboard">
        <table id="pipeline-table">
            <thead>
                <tr>
                    <th onclick="sortTable(0)">RCSB ID</th>
                    <th onclick="sortTable(1)">Status</th>
                    <th onclick="sortTable(2)">Duration</th>
"""

    # Add stage headers
    col_index = 3
    for stage in stages:
        stage_name = stage.replace('_', ' ').title()
        html += f"""                    <th onclick="sortTable({col_index})" class="status-cell">{stage_name}</th>
"""
        col_index += 1
    
    html += """                </tr>
            </thead>
            <tbody>
"""
    
    # Add rows for each structure
    for log in logs:
        rcsb_id = log.get('rcsb_id', 'Unknown')
        status = log.get('overall_status', 'pending')
        status_class = "success" if status == "success" else "failure" if status == "failure" else "pending"
        duration = log.get('duration', 0)
        
        html += f"""                <tr onclick="toggleRowSelection(this)" data-id="{rcsb_id}">
                    <td><strong>{rcsb_id}</strong></td>
                    <td class="{status_class}">{status.upper()}</td>
                    <td>{format_duration(duration)}</td>
"""
        
        # Add cells for each stage
        for stage in stages:
            stage_data = log.get('stages', {}).get(stage, {})
            stage_status = stage_data.get('status', 'pending')
            stage_class = "success" if stage_status == "success" else "failure" if stage_status == "failure" else "skipped" if stage_status == "skipped" else "pending"
            icon = "✓" if stage_status == "success" else "✗" if stage_status == "failure" else "⚠" if stage_status == "skipped" else "⋯"
            
            stage_duration = stage_data.get('duration', None)
            duration_str = format_duration(stage_duration) if stage_duration is not None else ""
            
            html += f"""                    <td class="status-cell {stage_class}" data-status="{stage_status}">
                        <div>
                            <span class="status-icon">{icon}</span>
                            <div class="duration">{duration_str}</div>
                        </div>
                    </td>
"""
        
        html += """                </tr>
"""
    
    html += """            </tbody>
        </table>
    </div>
    
    <!-- Selection panel -->
    <div id="selection-panel">
        <div id="selection-header">
            <div>
                <strong>Selected structures</strong>
                <span id="selection-count">0</span>
            </div>
            <div>
                <button onclick="copySelectedIds()" class="btn btn-primary">Copy to clipboard</button>
                <button onclick="clearSelectedRows()" class="btn btn-secondary">Clear</button>
            </div>
        </div>
        <div id="selection-text"></div>
    </div>

    <div class="footer">
        <p>Generated on """ + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + """ · Displaying """ + str(len(logs)) + """ structures</p>
    </div>
    
    <script>
        // ===== COLUMN SORTING =====
        var sortColumn = 0;
        var sortDirection = 1; // 1 for ascending, -1 for descending
        
        function sortTable(column) {
            const table = document.getElementById('pipeline-table');
            const headers = table.getElementsByTagName('th');
            const tbody = table.tBodies[0];
            const rows = Array.from(tbody.getElementsByTagName('tr'));
            
            // Update sort direction
            if (sortColumn === column) {
                sortDirection = -sortDirection;
            } else {
                sortDirection = 1;
                sortColumn = column;
            }
            
            // Remove sorting indicators from all headers
            for (let i = 0; i < headers.length; i++) {
                headers[i].classList.remove('sorted-asc', 'sorted-desc');
            }
            
            // Add sorting indicator to current header
            if (sortDirection === 1) {
                headers[column].classList.add('sorted-asc');
            } else {
                headers[column].classList.add('sorted-desc');
            }
            
            // Sort rows
            rows.sort((a, b) => {
                let aValue, bValue;
                
                // Handle special sorting for different column types
                if (column === 0) {
                    // Sort by RCSB ID
                    aValue = a.cells[column].textContent.trim();
                    bValue = b.cells[column].textContent.trim();
                } else if (column === 1) {
                    // Sort by status
                    aValue = a.cells[column].textContent.trim();
                    bValue = b.cells[column].textContent.trim();
                    
                    // Custom status ordering: SUCCESS, FAILURE, PENDING
                    const statusOrder = { 'SUCCESS': 0, 'FAILURE': 1, 'PENDING': 2 };
                    aValue = statusOrder[aValue] ?? 3;
                    bValue = statusOrder[bValue] ?? 3;
                } else if (column >= 3) {
                    // Sort by stage status
                    const aCell = a.cells[column];
                    const bCell = b.cells[column];
                    
                    const aStatus = aCell.getAttribute('data-status');
                    const bStatus = bCell.getAttribute('data-status');
                    
                    // Custom status ordering for stages: success, failure, skipped, pending
                    const stageOrder = { 'success': 0, 'failure': 1, 'skipped': 2, 'pending': 3 };
                    aValue = stageOrder[aStatus] ?? 4;
                    bValue = stageOrder[bStatus] ?? 4;
                } else {
                    // Default string comparison
                    aValue = a.cells[column].textContent.trim();
                    bValue = b.cells[column].textContent.trim();
                }
                
                // Compare and apply sort direction
                if (aValue < bValue) {
                    return -1 * sortDirection;
                } else if (aValue > bValue) {
                    return 1 * sortDirection;
                }
                return 0;
            });
            
            // Reorder rows
            rows.forEach(row => {
                tbody.appendChild(row);
            });
        }
        
        // ===== ROW SELECTION =====
        const selectedRows = new Set();
        const selectionPanel = document.getElementById('selection-panel');
        const selectionText = document.getElementById('selection-text');
        const selectionCount = document.getElementById('selection-count');
        
        function toggleRowSelection(row) {
            const id = row.getAttribute('data-id');
            
            if (row.classList.contains('selected')) {
                row.classList.remove('selected');
                selectedRows.delete(id);
            } else {
                row.classList.add('selected');
                selectedRows.add(id);
            }
            
            updateSelectionPanel();
        }
        
        function updateSelectionPanel() {
            if (selectedRows.size > 0) {
                selectionPanel.style.display = 'block';
                selectionCount.textContent = selectedRows.size;
                selectionText.textContent = Array.from(selectedRows).join('\\n');
            } else {
                selectionPanel.style.display = 'none';
            }
        }
        
        function clearSelectedRows() {
            document.querySelectorAll('tr.selected').forEach(row => {
                row.classList.remove('selected');
            });
            selectedRows.clear();
            updateSelectionPanel();
        }
        
        function copySelectedIds() {
            if (selectedRows.size === 0) return;
            
            const text = Array.from(selectedRows).join('\\n');
            navigator.clipboard.writeText(text)
                .then(() => {
                    alert(`Copied ${selectedRows.size} structure IDs to clipboard`);
                })
                .catch(err => {
                    console.error('Failed to copy: ', err);
                    alert('Failed to copy. Please select the text manually and copy.');
                });
        }
        
        // ===== SEARCH FILTERING =====
        document.getElementById('search-input').addEventListener('input', function() {
            const searchText = this.value.toLowerCase();
            const rows = document.querySelectorAll('#pipeline-table tbody tr');
            
            rows.forEach(row => {
                const id = row.getAttribute('data-id').toLowerCase();
                if (id.includes(searchText)) {
                    row.style.display = '';
                } else {
                    row.style.display = 'none';
                }
            });
        });
        
        // Initially sort by the first column (RCSB ID)
        sortTable(0);
    </script>
</body>
</html>
"""
    
    # Write HTML to file
    with open(output_path, 'w') as f:
        f.write(html)
    
    print(f"Dashboard generated at {output_path}")

def main():
    """Main function to parse arguments and generate dashboard."""
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} /path/to/logs/directory [output.html]")
        sys.exit(1)
    
    log_dir = Path(sys.argv[1])
    if not log_dir.exists() or not log_dir.is_dir():
        print(f"Error: {log_dir} is not a valid directory")
        sys.exit(1)
    
    output_file = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("dashboard.html")
    
    # Scan log directory
    log_data = scan_logs_directory(log_dir)
    
    if not log_data:
        print(f"No valid log files found in {log_dir}")
        sys.exit(1)
    
    # Generate HTML
    generate_html(log_data, output_file)

if __name__ == "__main__":
    main()

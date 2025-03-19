#!/usr/bin/python
import os
import sys
import json
from glob import glob
from pathlib import Path
from datetime import datetime

import jinja2

def load_json(file_path):
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        return {"error": str(e)}

def generate_report(json_files, sample_name):
    current_date = datetime.now().strftime("%Y-%m-%d")
    print(json_files)
    
    sections = {}
    for json_file in json_files:
        section_name = os.path.splitext(json_file)[0]
        sections[section_name] = load_json(json_file)["interpretation"]
    
    template = jinja2.Template("""
    <html>
    <head>
        <title>Report for {{ sample_name }}</title>
    </head>
    <body>
        <h1>Sample Report: {{ sample_name }}</h1>
        <h2>Date: {{ current_date }}</h2>
        {% for section, data in sections.items() %}
        <h3>Section: {{ section }}</h3>
        <div>{{ data }}</div>
        {% endfor %}
    </body>
    </html>
    """)
    
    return template.render(sample_name=sample_name, current_date=current_date, sections=sections)

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <json_file(s)>")
        sys.exit(1)

    sample = sys.argv[1]
    json_files = sys.argv[2:]
    
    report_html = generate_report(json_files, sample)
    
    with open(f"{sample}.html", "w") as f:
        f.write(report_html)

if __name__ == "__main__":
    main()

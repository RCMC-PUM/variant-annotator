#!/usr/bin/python

import sys
import json
from glob import glob
from os.path import join
from pathlib import Path
from collections import defaultdict
from datetime import datetime

import vcfpy
from jinja2 import Environment, FileSystemLoader


def clean_value(val):
    val = str(val)

    for char in ["[", "]", "'"]:
        val = val.replace(char, "")

    val = val.replace("|", " | ")
    return val


def clean_name(name):
    return name.strip().replace("_", "")


def constrain(text, n = 25):
    if len(text) > n:
        text = text[:n]
        text += "..."
    return text


def load_data(path):
    reader = vcfpy.Reader.from_path(path)
    parsed = defaultdict(dict)

    for cnt, record in enumerate(reader):
        info = {a: clean_value(b) for a, b in record.INFO.items()}
        
        pos = f"{record.CHROM}:{record.POS}"
        ref = constrain(record.REF)
        alt = " | ".join([f"{constrain(alt.value)} [{alt.type}]" for alt in record.ALT])

        gt = record.calls[0].__dict__["data"].get("GT", "")
        rs = " | ".join([f"rs{rs}" for rs in record.INFO.get("RS", [])])
        
        parsed[cnt].update({"VARIANT": {"POS": pos, 
                                        "REF": ref, 
                                        "ALT": alt, 
                                        "GT": gt, 
                                        "RS ID": rs}
                           })
        
        parsed[cnt].update({"INFO": info})
        parsed[cnt].update({"Calls": record.calls[0].__dict__["data"]})

    return parsed


def load_json(path):
    with open(path, "r") as handle:
        file = json.load(handle)
    return file


def render(data, sample, template):
    env = Environment(loader=FileSystemLoader("."))
    template = env.get_template(Path(template).name)

    # Render the template with data
    output_html = template.render(data)

    # Save the rendered HTML to a file
    with open(f"{sample}.html", "w") as f:
        f.write(output_html)

    print("HTML file generated successfully! Open 'output.html' to view it.")


def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <vcf_file(s)>")
        sys.exit(1)

    sample = sys.argv[1]

    callers = sys.argv[2].split(",")
    files = sys.argv[3].split(",")
    clinical_annotation =  sys.argv[4]

    workflow_params = load_json(sys.argv[5])
    annots_params = load_json(sys.argv[6])

    template = sys.argv[7]
    
    report_data = {"sample_name": sample,
                   "date": datetime.today().strftime('%Y-%m-%d'),
                   "workflow": workflow_params,
                   "annots": annots_params,
                   "clinical_annots": clinical_annotation,
                   "variants": {}, 
                  }
    
    zipped = dict(zip(callers, files))
    zipped_sorted = {k: zipped[k] for k in ["small_variant", "cnv", "repeats", "sv", "ploidy"]}
    
    for caller, file in zipped_sorted.items():
        data = load_data(file)
        report_data["variants"][caller] = data

    with open(f"{sample}.json", "w") as handle:
        json.dump(report_data, handle, indent=2)

    render(report_data, sample, template)


if __name__ == "__main__":
    main()

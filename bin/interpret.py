#!/usr/bin/python

from pathlib import Path
import time
import json
import sys

import pandas as pd
from openai import OpenAI

def main():
    # Check if the right number of arguments is provided
    if len(sys.argv) != 5:
        print("Usage: interpret.py <path_to_tsv_file> <caller_name> <clinical_data> <model_config>")
        sys.exit(1)

    # Get command line arguments
    tsv_file_path = sys.argv[1]
    caller_name = sys.argv[2]
    clinical_data = sys.argv[3]
    
    if not clinical_data.strip():
        clinical_data = "Not provided" 

    model_config = sys.argv[4]

     # Load the prompt from the prompt.json file
    with open(model_config, "r") as prompt_file:
        model_config = json.load(prompt_file)
        
        model = model_config.get("model", "")
        prompt = model_config.get("prompt", "")
        api_key = model_config.get("api_key", "")
        instructions = model_config.get("instructions", "")
        
    client = OpenAI(api_key=api_key)

    try:
        # Load the TSV file
        df = pd.read_csv(tsv_file_path, sep="\t")
        if not df.empty:
            # Prepare variant descriptions for analysis
            variants = " ".join([f"Variant {x};" for x in df.to_dict(orient="records")])
            # Format the prompt with the clinical data and variants
            formatted_prompt = f"{prompt} (1) clical data: {clinical_data}; (2) caller type: {caller_name}; (3) detected variants {json.dumps(variants, indent=2)}."
    
            # Send the prompt to OpenAI GPT
            response = client.responses.create(
                model=model,
                instructions=instructions,
                input=formatted_prompt,
            )
            
            # Interpretation result
            interpretation = response.output_text
            output_data = {
                    "interpretation": interpretation,
                    "file_name": tsv_file_path,
                    "caller_name": caller_name
                }
            time.sleep(1)
            
        else:
            output_data = {
                "interpretation": "No variants detected",
                "file_name": tsv_file_path,
                "caller_name": caller_name
            }

    except pd.errors.EmptyDataError:
        output_data = {
                "interpretation": "No variants detected",
                "file_name": tsv_file_path,
                "caller_name": caller_name
            }

    output_file = Path(tsv_file_path).name
    with open(f"{output_file}.json", "w") as handle:
        json.dump(output_data, handle, indent=2)

if __name__ == "__main__":
    main()
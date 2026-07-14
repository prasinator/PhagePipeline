import json
import os

try:
    with open('tools/zenodo.json', 'r') as f:
        data = json.load(f)
    print("Files found in Zenodo record:")
    for file_info in data.get('files', []):
        name = file_info.get('key')
        link = file_info.get('links', {}).get('self')
        size = file_info.get('size')
        print(f"- Name: {name}")
        print(f"  Link: {link}")
        print(f"  Size: {size} bytes")
except Exception as e:
    print(f"Error parsing JSON: {e}")

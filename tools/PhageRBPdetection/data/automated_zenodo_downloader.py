import urllib.request
import json
import os
import sys

# 1. Fetch live metadata to make sure we have the exact right download URL
url = "https://zenodo.org/api/records/10515367"
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})

print("Fetching metadata from Zenodo record 10515367...")
try:
    with urllib.request.urlopen(req) as response:
        data = json.loads(response.read().decode())
    
    # 2. Look for the ESM fine-tuned weights file
    target_file = "RBPdetect_v3_ESMfineT33.zip"
    download_url = None
    file_size = 0
    
    for f in data.get('files', []):
        if f.get('key') == target_file:
            # Zenodo stores the direct download link under 'content' or 'download' inside links
            download_url = f.get('links', {}).get('content') or f.get('links', {}).get('download')
            file_size = f.get('size', 0)
            break
            
    if not download_url:
        print(f"Error: Could not find {target_file} in this Zenodo record.")
        sys.exit(1)
        
    print(f"Found target file! Size: {file_size / (1024*1024):.2f} MB")
    print(f"Downloading from: {download_url}")
    
    # 3. Stream download with fake User-Agent to bypass the 403 Forbidden error
    download_req = urllib.request.Request(download_url, headers={'User-Agent': 'Mozilla/5.0'})
    
    with urllib.request.urlopen(download_req) as response, open(target_file, 'wb') as out_file:
        downloaded = 0
        block_size = 65536  # 64kb chunks for faster streaming
        
        while True:
            buffer = response.read(block_size)
            if not buffer:
                break
            downloaded += len(buffer)
            out_file.write(buffer)
            
            # Progress tracker
            if file_size > 0:
                percent = downloaded * 100 / file_size
                sys.stdout.write(f"\rDownload Progress: {percent:.2f}% ({downloaded / (1024*1024):.2f} MB)")
                sys.stdout.flush()
                
    print("\n\nSuccess! Verifying zip integrity...")
    # Quickly check the file structure
    stream = os.popen(f'unzip -t {target_file}')
    output = stream.read()
    print(output)

except Exception as e:
    print(f"\nAn error occurred during pipeline execution: {e}")

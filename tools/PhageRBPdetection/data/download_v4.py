import urllib.request
import os
import zipfile
import sys

url = "https://zenodo.org/api/records/14810759/files/RBPdetect_v4_ESMfine.zip/content"
output = "tools/PhageRBPdetection/data/RBPdetect_v4_ESMfine.zip"
dest_dir = "tools/PhageRBPdetection/data"

print(f"Downloading {url} to {output}...")
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})

try:
    with urllib.request.urlopen(req) as response, open(output, 'wb') as out_file:
        total_size = int(response.getheader('Content-Length', 0))
        downloaded = 0
        block_size = 65536
        while True:
            buffer = response.read(block_size)
            if not buffer:
                break
            downloaded += len(buffer)
            out_file.write(buffer)
            if total_size > 0:
                percent = downloaded * 100 / total_size
                sys.stdout.write(f"\rProgress: {percent:.2f}% ({downloaded / (1024*1024):.2f} MB)")
                sys.stdout.flush()
    print("\nDownload complete! Unzipping...")
    with zipfile.ZipFile(output, 'r') as zip_ref:
        zip_ref.extractall(dest_dir)
    print("Unzipping complete!")
    os.remove(output)
    print("Temporary zip removed.")
except Exception as e:
    print(f"\nAn error occurred: {e}")

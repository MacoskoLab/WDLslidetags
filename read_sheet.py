import sys
import gspread
import re
import subprocess
import csv

empty = ['X','',' ']
def checkgsfile(path):
    res = subprocess.run(['gsutil', 'ls', path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return True if res.returncode == 0 else False

#### Load Data #################################################################

# Get the input BCL
if len(sys.argv) < 2:
    print("Too few arguments, please input the BCL path")
    sys.exit()
if len(sys.argv) > 2:
    print("Too many argments, please input only the BCL path")
    sys.exit()
bcl = sys.argv[1]
print(f"BCL: {bcl}")
assert checkgsfile(f"gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/01_INBCLS/{bcl}")

# Get the sheet
gspread_client = gspread.service_account(filename="upload_for_google_key.json")
spreadsheet = gspread_client.open("Slide-Tags Experiment Log")
sheet = spreadsheet.worksheet("Tags")

# Get the data row matching the BCL
bclcol = sheet.find("BCL").col
rows = [cell.row for cell in sheet.findall(bcl) if cell.col == bclcol]
assert len(rows) <= 1, "input BCL had multiple spreadsheet matches"
assert len(rows) >= 1, "input BCL had no spreadsheet matches"
row = rows[0]

# Get the RNA Indexes
RNAindexes = sheet.cell(row, sheet.find("RNA Index").col).value.split('\n')
print(f"RNA Indexes: {RNAindexes}")

# Get the SB Indexes
SBindexes = sheet.cell(row, sheet.find("SB Index").col).value.split('\n')
print(f"SB Indexes: {SBindexes}")

# Get the Transcriptomes
transcriptomes = sheet.cell(row, sheet.find("Transcriptome").col).value.split('\n')
print(f"Transcriptomes: {transcriptomes}")
# If only one transcriptome, assume it applies to all RNA samples
if len(transcriptomes) == 1:
    transcriptomes *= len(RNAindexes)
print(f"Transcriptomes: {transcriptomes}")
# Check that the transcriptomes exist in the bucket
assert all(checkgsfile(f'gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/03_REFS/{ref}') for ref in transcriptomes if ref not in empty)

assert len(RNAindexes) == len(SBindexes) == len(transcriptomes) > 0
assert all(pair[1] not in empty for pair in zip(RNAindexes, transcriptomes) if pair[0] not in empty)

# Get the number of lanes
lanes = int(sheet.cell(row, sheet.find("Lanes").col).value)
bash = f"gsutil ls gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/01_INBCLS/{bcl}/Data/Intensities/BaseCalls | egrep '/L[0-9][0-9][0-9]/$'"
bcllanes = subprocess.check_output(bash, shell=True).decode('utf8').strip().split('\n')
assert lanes == len(bcllanes) > 0
print(f"Lanes: {lanes} identified")

# Get the RNA technique
RNAtech = sheet.cell(row, sheet.find("RNA Technique").col).value
print(f"RNA Technique: {RNAtech}")

# Get the demultiplex staus
demultiplex = sheet.cell(row, sheet.find("Demultiplex").col).value
print(f"Demultiplex: {demultiplex}")

# Get the pucks
pucks = sheet.cell(row, sheet.find("Puck ID").col).value.split('\n')
print(f"Pucks: {pucks}")
pucks = [puck if (puck[-4:]==".csv" or puck in empty) else f"{puck}.csv" for puck in pucks]
print(f"Pucks: {pucks}")
# Check that the pucks exist in the bucket
assert all(checkgsfile(f'gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/05_PUCKS/{puck}') for puck in pucks if puck not in empty)

if demultiplex == "NO":
    assert len(SBindexes) == len(pucks) > 0
    assert all((pair[0] in empty)==(pair[1] in empty) for pair in zip(SBindexes,pucks))

# bcl, RNAindexes, SBindexes, transcriptomes, lanes, RNAtech, demultiplex, pucks

#### Create Files ##############################################################

# Write the Indexes.csv file
print("Writing Indexes.csv...")
with open("Indexes.csv", mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Lane', 'Sample', 'Index'])
    [[writer.writerow([i+1, "RNA-"+ind, ind]) for i in range(lanes)] for ind in RNAindexes if ind not in empty]
    [[writer.writerow([i+1, "SB-"+ind, ind]) for i in range(lanes)] for ind in SBindexes if ind not in empty]
print("Done")

# Write the Transcriptomes.csv file
print("Writing Transcriptomes.csv...")
if RNAtech == "count":
    with open("Transcriptomes.csv", mode='w', newline='') as file:
        writer = csv.writer(file)
        [writer.writerow([RNAtech, pair[0], pair[1]]) for pair in zip(RNAindexes,transcriptomes) if pair[0] not in empty]
elif RNAtech == "multi (FFPE)":
    print(f"{RNAtech}: not yet implemented  - nothing to write")
else:
    print(f"{RNAtech}: not yet implemented  - nothing to write")
print("Done")

# Write the Spatial.csv file
print("Writing Spatial.csv...")
if demultiplex == "NO" and RNAtech == "count":
    m = [tup[0] not in empty and tup[1] not in empty for tup in zip(SBindexes, RNAindexes)]
    SBkeys      = [ind for ind,b in zip(SBindexes,m) if b]
    puckpaths   = [f"gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/05_PUCKS/{puck}" for puck,b in zip(pucks,m) if b]
    cbpaths     = [f"gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/04_COUNTS/{bcl}/RNA-{ind}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" for ind,b in zip(RNAindexes,m) if b]
    fullcbpaths = [f"gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/04_COUNTS/{bcl}/RNA-{ind}/outs/raw_feature_bc_matrix/barcodes.tsv.gz" for ind,b in zip(RNAindexes,m) if b]
    with open("Spatial.csv", mode='w', newline='') as file:
        writer = csv.writer(file)
        [writer.writerow([pair[i] for i in range(4)]) for pair in zip(SBkeys,puckpaths,cbpaths,fullcbpaths)]
else:
    print(f"Not yet implemented  - nothing to write")
print("Done")


import sys
import gspread
import re
import subprocess
import csv

if len(sys.argv) != 3:
    print("Usage: read_sheet.py BCL bucket")
    sys.exit()

print("-----RUNNING READ_SHEET.PY-----")

empty = ['X','',' ']
def checkgsfile(path):
    res = subprocess.run(['gsutil', 'ls', path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return True if res.returncode == 0 else False

#### Set DNS ###################################################################

# The internal DNS within a terra instance redirects googleapis.com -> restricted.google.com and bans the request
# Don't mess with /etc/resolv.conf, it seems to have side effects
# Modify the python to use dns_server=8.8.8.8 without affecting the rest (hopefully)
import dns.resolver
import socket
def custom_getaddrinfo(host, port, family=0, type=0, proto=0, flags=0, dns_server='8.8.8.8'):
    resolver = dns.resolver.Resolver()
    resolver.nameservers = [dns_server]
    try:
        answer_ipv4 = resolver.resolve(host, 'A')
    except (dns.resolver.NoAnswer, dns.resolver.NXDOMAIN, dns.resolver.Timeout):
        answer_ipv4 = []
    try:
        answer_ipv6 = resolver.resolve(host, 'AAAA')
    except (dns.resolver.NoAnswer, dns.resolver.NXDOMAIN, dns.resolver.Timeout):
        answer_ipv6 = []
    if not answer_ipv4 and not answer_ipv6:
        raise socket.gaierror(socket.EAI_NONAME, 'Name or service not known')
    addresses = [(socket.AF_INET, a.to_text()) for a in answer_ipv4] + [(socket.AF_INET6, a.to_text()) for a in answer_ipv6]
    results = []
    for af, addr in addresses:
        if family != 0 and family != af:
            continue
        if type == socket.SOCK_STREAM:
            proto = socket.IPPROTO_TCP
        elif type == socket.SOCK_DGRAM:
            proto = socket.IPPROTO_UDP
        else:
            raise ValueError("Unknown socket type")
        results.append((af, type, proto, '', (addr, port)))
    if not results:
        raise socket.gaierror(socket.EAI_FAMILY, 'Address family not supported')
    return results
# Replace the default getaddrinfo with custom_getaddrinfo
socket.getaddrinfo = custom_getaddrinfo

#### Open Sheet ################################################################

# Get the input BCL
bcl = sys.argv[1].replace("/", "").replace("\\", "") ; print(f"BCL: {bcl}")
bucket = sys.argv[2].replace("/", "").replace("\\", "") ; print(f"Bucket: {bucket}")
assert checkgsfile(f"gs://{bucket}"), f"ERROR: the path gs://{bucket} could not be accessed"

print(f"BCL exists:   {checkgsfile(f'gs://{bucket}/01_BCLS/{bcl}')}")
print(f"FASTQS exist: {checkgsfile(f'gs://{bucket}/02_FASTQS/{bcl}')}")
print(f"COUNTS exist: {checkgsfile(f'gs://{bucket}/03_COUNTS/{bcl}')}")

# Load the worksheets
gspread_client = gspread.service_account(filename="google_key.json")
spreadsheet = gspread_client.open("Slide-Tags Experiment Log")

# Search for the BCL in the spreadsheet
RNAtech = None
for tech in ["Tags", "FFPE"]:
    sheet = spreadsheet.worksheet(tech)
    rows = [cell.row for cell in sheet.findall(bcl) if cell.col == sheet.find("BCL").col]
    if len(rows) == 0: continue
    assert len(rows) <= 1, f"ERROR: input BCL {bcl} had {len(rows)} spreadsheet matches in {tech}"
    assert RNAtech == None, f"ERROR: input BCL {bcl} had multiple spreadsheet matches: found in {RNAtech}, {tech}"
    RNAtech = tech
assert RNAtech!=None, f"ERROR: input BCL {bcl} had no spreadsheet matches"
print(f"RNA Technique: {RNAtech}")

# Get the worksheet row matching the BCL
sheet = spreadsheet.worksheet(RNAtech)
row = [cell.row for cell in sheet.findall(bcl) if cell.col == sheet.find("BCL").col][0]
print(f"Row number: {row}")

#### Load Data #################################################################

# Get the RNA Indexes
RNAindexes = sheet.cell(row, sheet.find("RNA Index").col).value
if RNAindexes != None:
    RNAindexes = RNAindexes.split('\n')
    print(f"RNA Indexes: {RNAindexes}")
else:
    print("No RNA indexes found")

# Get the SB Indexes
SBindexes = sheet.cell(row, sheet.find("SB Index").col).value
if SBindexes != None:
    SBindexes = SBindexes.split('\n')
    print(f"SB Indexes: {SBindexes}")
else:
    print("No SB indexes found")

# Get the transcriptomes
transcriptomes = sheet.cell(row, sheet.find("Transcriptome").col).value
if transcriptomes != None:
    transcriptomes = transcriptomes.split('\n')
    if len(transcriptomes) == 1:
        transcriptomes *= len(RNAindexes)
    print(f"Transcriptomes: {transcriptomes}")
    assert len(transcriptomes) == len(RNAindexes) > 0, "Number of transcriptomes does not match the number of RNA indexes, add X for blank"
    assert all(checkgsfile(f'gs://{bucket}/references/{ref}') for ref in transcriptomes if ref not in empty), "Not all transcriptomes exist in the bucket"
else:
    print("No transcriptomes found")

# Get the probe sets
if RNAtech == "FFPE":
    probes = sheet.cell(row, sheet.find("Probe set").col).value
    if probes != None:
        probes = probes.split('\n')
        if len(probes) == 1:
            probes *= len(RNAindexes)
        print(f"Probe sets: {probes}")
        assert len(probes) == len(RNAindexes) > 0, "Number of probe sets does not match the number of RNA indexes"
        assert all(checkgsfile(f'gs://{bucket}/probesets/{ref}') for ref in probes if ref not in empty), "Not all probe sets exist in the bucket"
    else:
        print("No probe sets found")
else:
    probes = None

# Get the number of lanes - could skip spreadsheet and just pull from BCL if wanted
if checkgsfile(f"gs://{bucket}/01_BCLS/{bcl}"):
    lanes = sheet.cell(row, sheet.find("Lanes").col).value
    if lanes!=None and lanes.isdigit():
        lanes = int(lanes)
        bash = f"gsutil ls gs://{bucket}/01_BCLS/{bcl}/Data/Intensities/BaseCalls | egrep '/L[0-9][0-9][0-9]/$'"
        bcllanes = subprocess.check_output(bash, shell=True).decode('utf8').strip().split('\n') # if this step fails then no lanes in the BCL were found matching the pattern
        assert lanes == len(bcllanes), f"ERROR: spreadsheet lanes ({lanes}) did not match BCL lanes ({len(bcllanes)})"
        print(f"Lanes: {lanes} identified")
        lanes = None if lanes==0 else lanes
    else:
        print(f"No valid spreadsheet lanes found (input: {lanes})")
        lanes = None
else:
    print("No BCL found, skipping lane calculation")
    lanes = None

# Get the pucks as a list of lists
pucks = sheet.cell(row, sheet.find("Puck ID").col).value
if pucks!=None:
    pucks = [[puck.strip() for puck in pucklist.split(",")] for pucklist in pucks.split('\n')]
    pucks = [[f"{puck}.csv" if (puck[-4:]!=".csv") else puck for puck in pucklist if puck not in empty] for pucklist in pucks]
    assert len(pucks) == len(SBindexes), f"Number of pucks ({len(pucks)}) does not match the number of SB indexes ({len(SBindexes)}), add X for blank"
    for puck in [puck for pucklist in pucks for puck in pucklist]:
        if not checkgsfile(f'gs://{bucket}/pucks/{puck}'):
            assert False, f"Puck {puck} does not exist in the bucket, aborting..."
    print(f"Pucks: {pucks}")
else:
    pucks = []
    print("No pucks found")

# Get the demultiplex status
demultiplex = sheet.cell(row, sheet.find("Demultiplex").col).value
assert demultiplex in ["YES","NO"]
if demultiplex == "NO":
    assert len(RNAindexes) == len(SBindexes), "Demultiplex==NO, but the number of SBindexes and RNAindexes differs"
print(f"Demultiplex: {demultiplex}")

# Assert no commas, tabs, or spaces in any of the strings
allvals = [RNAindexes, SBindexes, transcriptomes, probes, *pucks]
allvals = [l if isinstance(l,list) else [l] for l in allvals]
allvals = [item for sublist in allvals for item in sublist]
assert all(val==None or (',' not in val) for val in allvals), "ERROR: Remove all commas from input"
assert all(val==None or ('\t' not in val) for val in allvals), "ERROR: Remove all tabs from input"
assert all(val==None or (' ' not in val) for val in allvals), "ERROR: Remove all spaces from input"

# Refuse to run mkfastq if mkfastq output already exists
if checkgsfile(f"gs://{bucket}/02_FASTQS/{bcl}"):
    print(f"mkfastq output for {bcl} already exists, will not run mkfastq")
    lanes = None

#### Create Files ##############################################################

run_mkfastq = (isinstance(RNAindexes,list) or isinstance(SBindexes,list)) and isinstance(lanes,int)
run_RNAcounts = isinstance(RNAindexes,list) and isinstance(transcriptomes,list) and (isinstance(probes,list) if RNAtech=="FFPE" else True)
run_SBcounts = isinstance(SBindexes,list) and isinstance(pucks, list)
run_make_seurat = demultiplex=="NO" and isinstance(RNAindexes, list) and isinstance(SBindexes, list)

# Write the Indexes.csv file
if run_mkfastq:
    print("Writing Indexes.csv...")
    with open("Indexes.csv", mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Lane', 'Sample', 'Index'])
        if isinstance(RNAindexes,list):
            [[writer.writerow([i+1, ind, ind]) for i in range(lanes)] for ind in RNAindexes if ind not in empty]
        if isinstance(SBindexes,list):
            [[writer.writerow([i+1, ind, ind]) for i in range(lanes)] for ind in SBindexes if ind not in empty]
    print(f"Done ({sum(1 for line in open('Indexes.csv'))} lines written)")
else:
    print(f"Not enough information to create Indexes.csv, writing a blank file...")
    open("Indexes.csv", 'w').close()

# Write the RNAcounts.tsv file
if run_RNAcounts:
    print("Writing RNAcounts.tsv...")
    probes = probes if probes!=None else ["X"]*len(RNAindexes)
    with open("RNAcounts.tsv", mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        for ind,ref,probe in zip(RNAindexes,transcriptomes,probes):
            if ind not in empty and ref not in empty:
                if not checkgsfile(f"gs://{bucket}/03_COUNTS/{bcl}/{ind}"):
                    if probe in empty:
                        writer.writerow([ind, ref])
                    else:
                        writer.writerow([ind, ref, probe])
                else:
                    print(f"{ind} already exists, skipping RNAcounts for this index")
    print(f"Done ({sum(1 for line in open('RNAcounts.tsv'))} lines written)")
else:
    print(f"Not enough information to create RNAcounts.tsv, writing a blank file...")
    open("RNAcounts.tsv", 'w').close()

# Write the SBcounts.tsv file
if run_SBcounts:
    print("Writing SBcounts.tsv...")
    with open("SBcounts.tsv", mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        for ind,pucklist in zip(SBindexes,pucks):
            if ind not in empty and len(pucklist)>0:
                writer.writerow([ind, ",".join(pucklist)])
    print(f"Done ({sum(1 for line in open('SBcounts.tsv'))} lines written)")
else:
    print(f"Not enough information to create SBcounts.tsv, writing a blank file...")
    open("SBcounts.tsv", 'w').close()

# Write the Spatial.tsv
if run_make_seurat:
    print("Writing Spatial.tsv...")
    with open("Spatial.tsv", mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        for RNAind,SBind in zip(RNAindexes, SBindexes):
            if RNAind not in empty and SBind not in empty:
                writer.writerow([f"gs://{bucket}/03_COUNTS/{bcl}/{RNAind}", f"gs://{bucket}/04_SPATIAL/{bcl}/{SBind}"])
    print(f"Done ({sum(1 for line in open('Spatial.tsv'))} lines written)")
else:
    print(f"Not enough information to create Spatial.tsv, writing a blank file...")
    open("Spatial.tsv", 'w').close()

print("-----COMPLETED READ_SHEET.PY-----")


### DOCUMENTATION ###
# This script reads the spreadsheet, validates the input, and writes files needed for downstream tasks
# See the workflow README for information on the input/output for these tasks

### ASSERTIONS ###
# The bucket must exist and be accessible
# 'upload_for_google_key.json' is required to access the spreadsheet
# The BCL must appear uniquely in the BCL column across all spreadsheets
# Number of transcriptomes must match the number of RNA indexes (add X for blank rows)
#     exception: if only one transcriptome is provided, it is assumed to apply to all indexes
#     also, all transcriptomes must exist in the bucket
# Same thing for the probes if doing FFPE
# Number of lanes in the sheet must match the number of lanes in the BCL Data/Intensities/BaseCalls (if provided)
# If pucks are provided, they must exist in the bucket and have the same length as SBindexes
# Demultiplex must be "YES" or "NO"
#   if NO, then #RNAindexes == #SBindexes and each has one puck (empty values allowed)
# No commas, tabs, or spaces allowed in RNAindexes, SBindexes, transcriptomes, probes, or pucks
#   exception is that commas are allowed in pucks

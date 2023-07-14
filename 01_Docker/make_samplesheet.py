import dns.resolver
import socket

# This is really really annoying
# But the internal DNS within a terra instance redirects googleapis.com ->
# restricted.google.com , and bans the request
# Tried a lot of different ways, but if mess with /etc/resolv.conf seems to have side effects
# This way, I just modify the python to use dns_server=8.8.8.8 without
# affecting the rest, hopefully
# ChatGPT helped with some encouragement
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


import sys
if len(sys.argv) < 2:
    print("Too few arguments, please input the BCL path")
    sys.exit()
if len(sys.argv) > 2:
    print("Too many argments, please input only the BCL path")
    sys.exit()
bcl = sys.argv[1]


import gspread
import re
gspread_client = gspread.service_account(filename="/upload_for_google_key.json")
sheets = gspread_client.open("Slide-Tags Experiment Log")
sheet = sheets.sheet1

rows = sheet.findall(re.compile('.*{}.*'.format(bcl)))
assert len(rows) <= 1, "input BCL had multiple spreadsheet matches"
assert len(rows) >= 1, "input BCL had no spreadsheet matches"
row = rows[0].row

col1 = sheet.find("RNA Index").col
col2 = sheet.find("Spatial Index").col

val1 = sheet.cell(row, col1).value.split('\n')
val2 = sheet.cell(row, col2).value.split('\n')
assert len(val1) == len(val2)
val1 = list(filter(lambda x: x != 'X' and x != "", val1))
val2 = list(filter(lambda x: x != 'X' and x != "", val2))

lanes = int(sheet.cell(row,sheet.find("Lanes").col).value)


import csv
with open("Indexes.csv", mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Lane', 'Sample', 'Index'])
    for val in val1:
        for i in range(lanes):
            writer.writerow([i+1, "RNA-"+val, val])
    for val in val2:
        for i in range(lanes):
            writer.writerow([i+1, "SB-"+val, val])

print("Done")

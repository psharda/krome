# read UMIST 2022 rates and convert them to KROME format

fname_umist = "rate22_final.rates"  # input file
fname = "network_umist.dat"  # output file


skip = ["PHOTON", "CRPHOT", "CRP"]
body = "@format:idx,R,R,P,P,P,P,tmin,tmax,rate\n"
body += "@common:user_av,user_crate\n"
count = 0
for row in open(fname_umist):
    srow = row.strip()
    if srow == "" or srow.startswith("#"):
        continue
    arow = srow.split(":")
    rtype = arow[1]
    rr = arow[2:4]
    pp = arow[4:8]
    ka, kb, kc = [float(x) for x in arow[9:12]]
    tmin, tmax = [float(x) for x in arow[12:14]]
    # print(rr, pp, ka, kb, kc)
    rate = None
    if rtype == "CR":
        rate = "%.2e * user_crate" % kc
    elif rtype == "PH":
        rate = "%.2e * exp(-%.2f * user_av)" % (ka, kc)
    else:
        rate = "%.2e" % ka
        if kb != 0e0:
            rate += " * (Tgas / 3e2)**(%.2f)" % kb
        if kc != 0e0:
            rate += " * exp(-%.2f / Tgas)" % kc

    body += "%d,%s,%s,%.2e,%.2e,%s\n" % (count, ",".join(rr), ",".join(pp), tmin, tmax, rate)
    count += 1

body = body.replace(",e-,", ",E,")

for s in skip:
    body = body.replace(",%s," % s, ",,")

with open(fname, "w") as f:
    f.write(body)

print("Wrote %d reactions to %s" % (count, fname))
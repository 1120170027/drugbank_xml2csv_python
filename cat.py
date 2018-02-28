from collections import defaultdict
for filename in ["./drug2carrier.csv", "./drug2enzyme.csv",
                 "./drug2target.csv", "./drug2transporter.csv"]:
    a = defaultdict(lambda: 0)
    ind = 1
    for line in open(filename):
        if ind != 1:
            for ty in line.strip().split(",")[2:]:
                a[ty] += 1
        else:
            pass
        ind += 1

    print("\n####### {} #######\n".format(filename))
    for cat, num in sorted(a.items(), key=lambda t: t[1], reverse=True):
        print cat, num

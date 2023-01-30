import glob

DIR = "./result/*"
files = glob.glob(DIR)

num = len(files)

d = {
    "INITIAL_COUNT" : 0,
    "LOOP_COUNT" : 0,
    "Score" : 0
}
for file in files:
    f = open(file, 'r')

    rows = f.readlines()

    for row in rows:
        a, b, c = row.split()
        d[a] += int(c)

print("INITIAL_COUNT  = " + str(d["INITIAL_COUNT"] / num))
print("LOOP_COUNT     = " + str(d["LOOP_COUNT"] / num))
print("TOTAL_SCORE    = " + str(d["Score"]))
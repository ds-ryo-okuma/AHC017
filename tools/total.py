import glob

DIR = "./result/*"
files = glob.glob(DIR)

count, total = 0, 0
for file in files:
    f = open(file, 'r')

    rows = f.readlines()

    a, b, c = rows[0].split()
    count += int(c)

    a, b, c = rows[1].split()
    total += int(c)

average = count / len(files)
print("LOOP_COUNT  = " + str(average))
print("TOTAL_SCORE = " + str(total))
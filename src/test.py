filename = "../input/to_download.txt"
with open(filename, "r") as f:
    next(f)
    for line in f.readlines():
        if line.rstrip() != "" and not line.startswith("#"):
            items = line.strip().split()
            if len(items) < 7:
            	print(line+"\n")
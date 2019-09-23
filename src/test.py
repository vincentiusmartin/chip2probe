import os
try:
	os.path.exists("../result/chipseq/")
	print([name for name in os.listdir(".") if os.path.isdir(name)])
except FileNotFoundError:
	print("Directory ../result/chipseq/ does not exist.") 
import argparse

from heterotypic_dictinput_example import inlist


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Run the pipeline.')
    parser.add_argument('index', type=int, help="Input index from dictinput.py")
    parser.add_argument('-o','--outdir', type=str, help="directory to output the results", default=".")
    args = parser.parse_args()

    findex = args.index
    outdir = "%s/%s_%s" % (args.outdir, inlist[findex]["tf1_name"], inlist[findex]["tf2_name"])
    

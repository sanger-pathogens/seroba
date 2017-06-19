import os

class Error (Exception): pass

def summarise(output_folder):
    folders = os.listdir(output_folder)
    with open('summary.tsv','w') as fobj:
        for d in folders:
            pred = os.path.join(output_folder,d,'pred.tsv')
            if os.path.isfile(pred):
                with open(pred, 'r') as f:
                    first_line = f.readline()
                fobj.write(first_line)

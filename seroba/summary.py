import os

class Error (Exception): pass

def summarise(output_folder):
    folders = os.listdir(folders)
    with open(os.path.join('summary.tsv'),'w') as fobj:
        for d in folders:
            pred = os.patj.join(output_folder,'pred.tsv')
            if os.path.isfile(pred):
                with open(pred, 'r') as f:
                    first_line = f.readline()
                fobj.write(first_line)

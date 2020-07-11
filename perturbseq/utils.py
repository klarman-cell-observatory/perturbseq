import sys

def display_progress(cur_num,max_num):
    sys.stdout.write('\r'+str(int(100*(1.0*cur_num/max_num)))+' %')
    sys.stdout.flush() 

def end_progress(num_spaces,newline=False):
    text="\r"
    for i in range(num_spaces):
        text=text+' '
    if newline:
        text=text+'\n'
    sys.stdout.write(text)

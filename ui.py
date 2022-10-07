import os

def clear_cmd_terminal():
    if os.name == "nt":
        os.system("cls")
    else:
        os.system("clear")

def generate_progress_bar(percent, text=True, width=20):
    ret_str = "["

    i = width
    for w in range(int((percent/100) * width)):
        ret_str += "#"
        i -= 1
        
    while i > 0:
        ret_str += "_"
        i -= 1

    ret_str += "]"

    if text:
        percent = round(percent, 2)
        ret_str += " " + str(percent) + "%"

    return ret_str

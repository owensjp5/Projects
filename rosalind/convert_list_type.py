def string_to_int_list(s):
    return [int(x) for x in s.split()]

def int_list_to_string(lst):
    return " ".join(str(x) for x in lst)
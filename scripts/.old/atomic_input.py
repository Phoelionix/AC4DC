import sys

command = sys.argv[1]

commands = {
'pulse': {},
'atom': {},
'numerical': {},
'output': {}
}

for cmd in commands:
    commands[cmd] += {}



def parse():
    try:
        commands[sys.argv[1]]
    except KeyError:
        

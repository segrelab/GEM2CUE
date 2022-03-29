# How to get the other item from the list

list = ['A', 'B']

choice = 'A'

other = [item for item in list if item != choice]

print(other)
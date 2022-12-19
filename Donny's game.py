import random

# first of all, all of your code should be under one function

def donny_game():
    first_number = random.randint(1,6)
    second_number = random.randint(6,12)

    # I prefer using f-strings for `Print` function

    print(f'First number is {first_number} and Second number is {second_number}')

    guess1 = int(input('Enter your number: '))
    
    # remember the input is in string format so you've to change it to integer data format.
    # a string answer =! integer answer even if the numbers are the same.

    correct_guess = first_number*second_number

    if guess1 == correct_guess:
        print('CORRECT')         #if you missed to convert the 'input' answer into Integer data
                                 # format, this answer will compare an integer and a string and 
                                 # will not give us a correct answer

    else:
        print('WRONG, Please retry')


donny_game()

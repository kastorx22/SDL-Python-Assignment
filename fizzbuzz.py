#Kateryna
import argparse
import sys


def fizzbuzz(limit, rules):
    #runs the fizzbuzz until it reaches the limit, using teh specified rules
    for i in range(1, limit+1):
        output = "".join(word for divisor, word in rules.items() if i % divisor == 0)
        print(output or i)


def limit_and_rules():
    parser=argparse.ArgumentParser(description='customizable fizzbuzz')
    #allow user to choose to the maximum value, default limit is 100
    parser.add_argument('--limit', type=int, default=100, help="the maximum number")
    #next 3 lines done by Ana, allows user to add custom rules and appends them to the existing set of rules
    parser.add_argument("--rule", action="append", nargs=2, metavar=("Factor", "Replacement Word"),
                        help="Add a factor and its replacement word")
    args=parser.parse_args()
    #dictionary of existing rules
    rules = {3: "Fizz", 5: "Buzz", 7: "Fang", 11: "Bang"}
    #next 3 lines done by Ana, updates the rules
    if args.rule:
        for factor, word in args.rule:
            rules[int(factor)]=word
    #Kateryna, returns the new dictionary of rules
    return args.limit, rules


if __name__ == "__main__":
    limit, rules = limit_and_rules()
    fizzbuzz(limit, rules)

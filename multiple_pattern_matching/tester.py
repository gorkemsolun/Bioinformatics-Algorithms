import argparse
import os
import re
import sys


def main():
    # Search for a pattern in a string which is taken from a file.
    pattern = "ACTAG"
    pattern_locations_to_check = []

    # open the file which contains the string to be searched
    with open(
        "C:\\Users\\gorke\\OneDrive\\CS_481\\Bioinformatics-Algorithms\\hw1\\short.fa",
        "r",
    ) as file:
        string = file.read()
        # get the second line
        string = string.split("\n")[1]

        # Using brute force check if the pattern is in the string
        # and print the index of the pattern in the string and control the index locations given
        matches = []
        for i in range(len(string)):
            is_match = string[i : i + len(pattern)] == pattern
            if is_match:
                matches.append(i)

        print("length of searched indexes: ", len(pattern_locations_to_check))
        print("length of found indexes: ", len(matches))

        for i in range(len(matches)):
            # print("Index of pattern in string: ", matches[i])
            if matches[i] != pattern_locations_to_check[i]:
                print("Wrong index found")
                print("Expected: ", pattern_locations_to_check[i])
                print("Found: ", matches[i])
                print("length of searched indexes: ", len(pattern_locations_to_check))
                print("length of found indexes: ", len(matches))


if __name__ == "__main__":
    main()

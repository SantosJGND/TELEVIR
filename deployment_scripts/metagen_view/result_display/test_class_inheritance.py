#!/usr/bin/python

from dataclasses import dataclass


class SimpleClass:
    def __init__(self, name: str, age: int):
        self.name = name
        self.age = age

    def __str__(self):
        print(f"name {self.name}")
        print(f"age {self.age}")


class Operations(SimpleClass):
    def double_age(self):
        return self.age * 2


if __name__ == "__main__":
    nobj = SimpleClass("steve", 20)
    print(SimpleClass.double_age())

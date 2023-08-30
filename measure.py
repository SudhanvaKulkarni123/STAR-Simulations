#this file contains some of the data structures that are 
#supposed to make the task of simulating physical systems easier



class Units:
    def __init__(self, name, power, magnitude):
        self.name = name
        self.power = power
        self.magnitude = 0
    def __str__(self):
        if self.power < 0:
            return str(self.magnitude) + "1/" + str(self.name) + "^" + str(self.power)
        if self.power > 0:
            return  str(self.magnitude) + str(self.name) + "^" + str(self.power)
        if self.power == 0:
            return 1
        


class Time(Units):
    def __init__(self, name, magnitude):
        if name == "sec" or name == "s" or name == "second" or name == "seconds":
            self.name = "s"
        if self.name == "min" or self.name == "minutes" or self.name == "minute":
            self.name = "min"
    def to_SI(self):
        if self.name == "min" :
            self.name = "s"
            self.magnitude = self.magnitude*60
        if self.name == "hr":
            self.name = "s"
            self.magnitude = self.magnitude*3600

        
class Length(Units):
    def __init__(self, name, magnitude):
        






        
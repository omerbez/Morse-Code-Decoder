from scipy.io import wavfile as wfl
from matplotlib import pyplot as plt
import os


class MorseCodeAnalyzer:

    morseLookupTable = {
        ".-": "A", "-...": "B", "-.-.": "C", "-..": "D", ".": "E", "..-.": "F", "--.": "G",
        "....": "H", "..": "I", ".---": "J", "-.-": "K", ".-..": "L", "--": "M", "-.": "N", "---": "O",
        ".--.": "P", "--.-": "Q", ".-.": "R", "...": "S", "-": "T", "..-": "U", "...-": "V", ".--": "W",
        "-..-": "X", "-.--": "Y", "--..": "Z", ".----": "1", "..---": "2", "...--": "3", "....-": "4",
        ".....": "5", "-....": "6", "--...": "7", "---..": "8", "----.": "9", "-----": "0",
        ".-.-.-": ".", "--..--": ",", "..--..": "?", "-..-.": "/", ".--.-.": "@", "---...": ":",
        "-.--.-": ")", ".-..-.": "\"", ".-.-.": "+", "-...-": "=", "-.--.": "(", "-....-": "-"
    }

    def __init__(self, path):

        if not os.path.isfile(path):
            raise FileExistsError("File doesn't exists!")

        (rate, sig) = wfl.read(path)

        self.__signal = (sig / 256.0)*2 + (-1)
        self.__NOISE_THRESHOLD = 10  # amount of sequential zeroes required to be considered as "no-signal"


    def __getSegmentsLengths(self):
        temp = []
        binarySignal = self.__signal.copy()

        # turn into binary signal, the signal may contain a little "noises" - some zero values
        # that should be 1. so we should handle that.
        binarySignal[binarySignal != 0] = 1
        i = 0
        while i < len(binarySignal):
            j = i+1
            while j < len(binarySignal):
                # splitted the while condition for more readability..
                if binarySignal[j] == binarySignal[i] or (self.__isANoise(binarySignal, binarySignal[i], j)):
                    j += 1
                else:
                    break

            if (j-i) > self.__NOISE_THRESHOLD:
                temp.append((binarySignal[i], j-i))
            i = j

        return temp


    def __isANoise(self, binarySignal, currentValue, index):
        """
        check if binarySignal[index] is just a noise or a new segment.
        new segment occurred when number if next elements are different than "currentValue".
        check ahead at most "NOISE_THRESHOLD" amount of elements.
        :param binarySignal: the converted signel
        :param currentValue: the current proceed value (0 or 1)
        :param index: the index we "encountered" in - binarySignal[index] != currentValue.
        :return: True if the current index is just a noise and it should be considered as "currentValue".
        """
        j = index + 1
        counter = 0
        while j < len(binarySignal) and binarySignal[j] != currentValue and counter < self.__NOISE_THRESHOLD:
            j += 1
            counter += 1

        return counter < self.__NOISE_THRESHOLD


    def __relationOf(self, n1, n2):
        if n1 < n2:
            return self.__relationOf(n2, n1)

        return n1 / float(n2)


    def __calcDotLength(self, a, b, rel):
        # rel can be 7/3, 7/1, 3/1

        if b > a:
            return self.__calcDotLength(b, a, rel)

        # if the relation is 3/1 or 7/1
        if (2.5 <= rel <= 3.9) or (6.0 <= rel <= 7.9):
            return b

        # relationship is 7/3.. dot size is a/7 or b/3..
        if 2.0 <= rel < 2.5:
            return a//7  # or b//3..

        # should not reach here..
        raise ValueError("Error: Couldn't calculate dot length. invalid rel argument")


    def __findDotLength(self, rawData):
        for i in range(len(rawData)):
            for j in range(i+1, len(rawData)):
                rel = self.__relationOf(rawData[i][1], rawData[j][1])
                if rel >= 2:
                    return self.__calcDotLength(rawData[i][1], rawData[j][1], rel)

        raise ValueError("Error: Couldn't calculate dot length. the signal might violate Morse-Code structure")


    def __calcMorseCodeText(self, rawData, dotLength):
        code = ""
        for seg in rawData:
            r = seg[1] / float(dotLength)
            if seg[0] == 1:
                if r <= 1.5:
                    code += "."
                elif 2.5 <= r <= 3.9:
                    code += "-"
                else:
                    raise ValueError("Error: Couldn't decode rawData into text. segment out of range! \n"
                                     "The signal might violate Morse code rules.")
            else:
                if r <= 1.5:
                    continue  # no-space - it's same letter..
                if 2.5 <= r <= 3.9:
                    code += " "  # space between letters
                elif 6.0 <= r <= 7.9:
                    code += "/"  # space between words.. (real space)

        return code


    def __decodeMorseMsg(self, code):
        i = 0
        msg = ""
        while i < len(code):
            if code[i] == " ":
                # next letter..
                i += 1
                continue
            elif code[i] == "/":
                # next word..
                msg += " "
                i += 1
                continue

            # else, look for the letter..
            j = i+1
            while j < len(code) and code[j] != " " and code[j] != "/":
                j += 1

            try:
                msg += MorseCodeAnalyzer.morseLookupTable[code[i:j]]
            except KeyError:
                # ignore.. it might be some special key which is not exists in our lookup table...
                pass

            i = j

        return msg


    def startDecode(self):
        rawData = self.__getSegmentsLengths()
        dotLength = self.__findDotLength(rawData)
        codeText = self.__calcMorseCodeText(rawData, dotLength)
        print(codeText)
        text = self.__decodeMorseMsg(codeText)
        print(text)


    def plotSignal(self):
        plt.figure("Morse Code")
        plt.plot(self.__signal)
        plt.show()


    @classmethod
    def getLookupTable(cls):
        return cls.morseLookupTable

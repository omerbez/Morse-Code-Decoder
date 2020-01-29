from morse_decoder import MorseCodeAnalyzer

analyzer = MorseCodeAnalyzer("data/Example 6.wav")
analyzer.startDecode()
analyzer.plotSignal()

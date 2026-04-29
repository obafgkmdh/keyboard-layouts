import requests, json

URL = "https://github.com/monkeytypegame/monkeytype/raw/refs/heads/master/frontend/static/quotes/english.json"

quotes = requests.get(URL).json()["quotes"]

with open("quotes.txt", "w") as f:
    for quote in quotes:
        f.write(quote["text"])
        f.write("\n")

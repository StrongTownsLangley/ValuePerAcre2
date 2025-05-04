import json

class JsonFileLoader:

    def load_json_file(file_path):
        """Load a JSON file with UTF-8 encoding to avoid character decoding issues."""
        with open(file_path, 'r', encoding='utf-8') as file:
            return json.load(file)
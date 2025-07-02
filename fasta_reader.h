#pragma once
#include <fstream>
#include <vector>
#include <string>

inline std::pair<std::vector<std::string>, std::vector<std::string>> readFastaWithNames(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::string> sequences;
    std::vector<std::string> names;

    std::string line, currentSeq, currentName;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        
        if (line.compare(0, 3, "\xEF\xBB\xBF") == 0)
            line = line.substr(3);

        if (line[0] == '>') {
            if (!currentSeq.empty()) {
                sequences.push_back(currentSeq);
                names.push_back(currentName);
                currentSeq.clear();
            }
            currentName = line.substr(1);  
        } else {
            for (char c : line)
                if (!isspace(c)) currentSeq += c;
        }
    }

    if (!currentSeq.empty()) {
        sequences.push_back(currentSeq);
        names.push_back(currentName);
    }

    return {sequences, names};
}


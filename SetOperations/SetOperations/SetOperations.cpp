#include <iostream>
#include <set>
#include "allsymbols.h"
#include "firstTextSymbols.h"
using namespace std;

int main()
{
    string a;
    string b;
    std::cout << "Enter first text:\n";
    getline(cin, a);
    std::cout << "Enter second text:\n";
    getline(cin, b);

    set<char> all = allSymbols(a, b);

    std::cout << "All symbols in the texts:\n";
    for (auto symbol : all)
        cout << symbol << " ";
    cout << "\n";


    set<char> firstTextOnly = firstTextSymbols(a, b);
    std::cout << "Symbols in the first text but not in the second:\n";
    for (auto s : firstTextOnly)
        cout << s << " ";
    cout << "\n";


    set<char> firstText = textSymbols(a);
    set<char> secondText = textSymbols(b);

    std::cout << "Symbols in first text:\n";
    for (auto s : firstText)
        cout << s << " ";
    cout << "\n";
    std::cout << "In first text " << firstText.size() << " symbols\n";

    std::cout << "Symbols in second text:\n";
    for (auto s : secondText)
        cout << s << " ";
    cout << "\n";
    std::cout << "In second text " << secondText.size() << " symbols\n";
}

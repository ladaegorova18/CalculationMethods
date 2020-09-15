#include "firstTextSymbols.h"

set<char> firstTextSymbols(string text1, string text2)
{
	set<char> symbols1 = textSymbols(text1);
	set<char> symbols2 = textSymbols(text2);
	set<char> symbols3;
	symbols3.insert(symbols1.begin(), symbols1.end());
	for (auto symbol : symbols2)
	{
		if (symbols3.find(symbol) != symbols3.end())
			symbols3.erase(symbol);
	}
	return symbols3;
}
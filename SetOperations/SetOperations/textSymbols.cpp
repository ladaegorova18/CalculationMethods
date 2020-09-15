#include "textSymbols.h"

set<char> textSymbols(string s)
{
	set<char> symbols;
	for (auto symbol : s)
	{
		if ((int)symbol >= 0 && (int)symbol <= 255)
			symbols.insert(symbol);
	}
	return symbols;
}
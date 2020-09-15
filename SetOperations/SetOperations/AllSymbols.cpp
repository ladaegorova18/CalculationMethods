#include "allsymbols.h"

set<char> allSymbols(string text1, string text2)
{
	set<char> symbols1 = textSymbols(text1);
	set<char> symbols2 = textSymbols(text2);
	set<char> symbols3;
	symbols3.insert(symbols1.begin(), symbols1.end());
	symbols3.insert(symbols2.begin(), symbols2.end());
	return symbols3;
}

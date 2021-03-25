

#include <string>
#include <iostream>

int main()
{
	std::string str1 = "To be or not to be, that is the question";
	std::string str2 = "only ";
	std::string str3 = str1.substr(6, 12);
	
	str1.insert(32, str2);
	str1.replace(str1.find("to be", 0), 5, "to jump");
	str1.erase(9, 4);
	
	std::cout << str1 << std::endl;
	  
	for(int i = 0; i < str3.length(); i++)
		std::cout << str3[i] << std::endl;
	
	return 0;
}


#How to show only next line after the matched one?
awk '/blah/{getline; print}' logfile
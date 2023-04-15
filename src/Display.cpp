/**
 * @file Display.cpp
 * @brief 
*/

/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/

#include <Display.h>
#include <csignal>
#include <iostream>

// Initialise static variables.
WINDOW* Display::win; std::string Display::header;

void Display::create_screen(){
    initscr();
    clear();
    noecho();    
    cbreak();
    int startx = (80 - WIDTH) / 2;
    int starty = (24 - HEIGHT) / 2;
    signal(SIGINT,Display::signalHandler);  // clean up for interrupt
    //win = newwin(HEIGHT, WIDTH, starty, startx);
    win = newwin(0, 0, 0, 0);
    box(win, 0 , 0);	
    wrefresh(win);
    keypad(win, TRUE);
    nodelay(win,TRUE); // don't wait for input  
}

/// Screen displays the contents of the stream, and only the contents.
void Display::show(const std::stringstream& spooky_stream){
    box(win, 0 , 0);
    waddstr(win,spooky_stream.str().c_str());
    wrefresh(win);   
    werase(Display::win);
}
void Display::show(const std::stringstream& spooky_stream,const std::stringstream& second_stream, bool do_erase){
    box(win, 0 , 0);
    waddstr(win,(spooky_stream.str()+second_stream.str()).c_str());
    wrefresh(win);   
    if (do_erase) werase(Display::win);
}
// void Display::deactivate(){
//     touchwin(stdscr);
// }
// void Display::reactivate(){
//     touchwin(win);
// }    
void Display::close(){
    clrtoeol();
    refresh();
    endwin();     
}  

// cleans up the terminal on interrupt
void Display::signalHandler( int signum ) {
    endwin();
    std::cout << "Window ended successfully after interrupt signal (" << signum << ") received.\n";
    std::exit(signum);  
}    

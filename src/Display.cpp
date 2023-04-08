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


WINDOW* Display::win; // The single line that demands a src file.

void Display::create_screen(){
    initscr();
    clear();
    noecho();    
    cbreak();
    int startx = (80 - WIDTH) / 2;
    int starty = (24 - HEIGHT) / 2;
    win = newwin(HEIGHT, WIDTH, starty, startx);
    refresh();
    keypad(win, TRUE);
    nodelay(win,TRUE); // don't wait for input      
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

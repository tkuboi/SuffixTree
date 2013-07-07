package bioGUI;

import java.io.*;
import java.util.ArrayList;
import java.lang.IndexOutOfBoundsException;

public class Position {
   public int idx;
   public int position;

   public Position(int i, int p) {
      idx = i;
      position = p;
   }

   public String toString() {
      return Integer.toString(idx) + ":" + Integer.toString(position);
   }
}

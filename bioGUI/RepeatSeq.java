package bioGUI;

import java.io.*;
import java.util.ArrayList;

public class RepeatSeq {
   public String sequence;
   public int count;
   public float expected;
   public ArrayList<Position> positions;

   public RepeatSeq() {
      this.sequence = "";
      this.count = 0;
      this.expected = 0;
      positions = new ArrayList<Position>();
   }

   public RepeatSeq( String seq, float expected, int count ) {
      this.sequence = seq;
      this.expected = expected;
      this.count = count;
      this.expected = 0;
      positions = new ArrayList<Position>();
   }

   public RepeatSeq( String seq, float expected, int count, ArrayList<Position> positions ) {
      this.sequence = seq;
      this.expected = expected;
      this.count = count;
      this.expected = 0;
      this.positions = positions;
   }

   public String toString() {
      String list = "";
      /*for (Position p : positions)
         list = list + p.toString() + ", ";
      if (list.length() > 3)
         list = list.substring(0, list.length() - 2);*/
      //return sequence + " : count= " + Integer.toString(count) + " Positions{" + list + "}";
      return String.format("%20s,          %5.2f, %8d\n",sequence,expected,count);
   }
}


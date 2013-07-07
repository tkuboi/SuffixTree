package bioGUI;

import java.io.*;
import java.util.ArrayList;
import java.lang.IndexOutOfBoundsException;

public class Node {
   private long id;
   private Position start;
   private Position end;
   //private int length;
   private int type;
   private Node[] children;
   private Node parent;
   private Node[] suffixLink;
   //private int position;
   private ArrayList<Position> positions;


   private static long currentId = 0;

   public static final int TYPE_BOTTOM = -1;
   public static final int TYPE_ROOT = 0;
   public static final int TYPE_INTERNAL = 1;
   public static final int TYPE_LEAF = 2;


   public Node() {
      id = currentId++;
      start = null;
      end = null;
      //length = 0;
      parent = null;
      type = TYPE_LEAF;
      this.children = new Node[6]; //place holder for A,T,C,G,$.
      this.suffixLink = new Node[6]; //place holder for A,T,C,G,$.
      this.positions = new ArrayList<Position>();
   }

   public Node(int type) {
      id = currentId++;
      start = null;
      end = null;
      //length = 0;
      parent = null;
      this.type = type;
      this.children = new Node[6]; //place holder for A,T,C,G,$.
      this.suffixLink = new Node[6]; //place holder for A,T,C,G,$.
      this.positions = new ArrayList<Position>();
   }

   public Node(int type, Position start, Position end) {
      id = currentId++;
      parent = null;
      this.start = start;
      this.end = end;
      //this.length = end.position - start.position + 1;
      this.type = type;
      this.children = new Node[6]; //place holder for A,T,C,G,$.
      this.suffixLink = new Node[6]; //place holder for A,T,C,G,$.
      this.positions = new ArrayList<Position>();
   }

   public Node(int type, Position start, Position end, Position pos) {
      id = currentId++;
      parent = null;
      this.start = start;
      this.end = end;
      //this.length = end.position - start.position + 1;
      this.type = type;
      this.children = new Node[6]; //place holder for A,T,C,G,$.
      this.suffixLink = new Node[6]; //place holder for A,T,C,G,$.
      this.positions = new ArrayList<Position>();
      this.positions.add(pos);
   }

   public long getId() {
      return this.id;
   }
   public void setStart(Position start) {
      this.start = start;
   }

   public Position getStart() {
      return this.start;
   }

   public void setEnd(Position end) {
      this.end = end;
   }

   public Position getEnd() {
      return this.end;
   }

   public void setType(int type) {
      this.type = type;
   }

   public int getType() {
      return this.type;
   }

   public void resetLabel(Position start, Position end) {
      this.start = start;
      this.end = end;
      //this.length = end.position - start.position + 1;
   }

   public void resetStart(Position start) {
      this.start = start;
      //this.length = end.position - start.position + 1;
   }

   public void resetStart(int n) {
      this.start.position += n;
      //this.length = end.position - start.position + 1;
   }

   public int getLength() {
      return end.position - start.position + 1;
   }

   public int getNumChildren() {
      int count = 0;
      for(int i = 0; i < 6; i++) {
         if(this.children[i] != null)
            count++;
      }
      return count;
   }

   /*public int addChild(Node child) {
      this.children.add(child);
      return this.children.size();
   }*/

   public void addChildAt(int idx, Node child) {
      this.children[idx] = child;
   }

   public Node getChildAt(int i) {
      //System.out.println(i);
      return this.children[i];
   }
    
   public Node removeChildAt(int i) {
      Node removed = null;
      removed = this.children[i];
      this.children[i] = null;
      return removed;
   }

   public ArrayList<Node> getChildren() {
      ArrayList<Node> list = new ArrayList<Node>();
      for(int i = 0; i < 6; i++){
         list.add(this.children[i]);
      }
      return list;
   }

   public int getNumPositions() {
      return this.positions.size();
   }

   public Position getPositionAt(int i) {
      return this.positions.get(i);
   }

   public int addPosition(Position pos) {
      this.positions.add(pos);
      return this.positions.size();
   }

   public int addPositions(ArrayList<Position> positions) {
      this.positions.addAll(positions);
      return this.positions.size();
   }

   public ArrayList<Position> getPositions() {
      return this.positions;
   }

   public void setParent(Node parent) {
      this.parent = parent;
   }

   public Node getParent() {
      return this.parent;
   }

   /*public void setSuffixLink(Node link) {
      this.suffixLink = link;
   }*/

   public void setSuffixLink(Node link, int i) {
      this.suffixLink[i] = link;
   }

   /*public Node getSuffixLink() {
      return this.suffixLink;
   }*/

   public Node getSuffixLink(int i) {
      return this.suffixLink[i];
   }

   public static void main(String[] args){
      Node node = new Node();
      System.out.println(node.getId());
      Node node2 = new Node();
      System.out.println(node2.getId());
   }
}

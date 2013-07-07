package bioGUI;

import java.io.*;
import java.util.ArrayList;
import java.io.IOException;
import java.lang.Math;

public class SuffixTree {
   //private String m_dna;
   private Node root;
   private Node leaf;
   private ArrayList<String> m_sequences;
   private ArrayList<Node> m_leaves;
   private float totalA;
   private float totalC;
   private float totalG;
   private float totalT;
   private int totalNucleotides;
   
   public SuffixTree() {
      root = null;
      m_sequences = new ArrayList<String>();
      m_leaves = new ArrayList<Node>();
   }

   public SuffixTree(String str) {
      root = null;
      m_sequences = new ArrayList<String>();
      m_sequences.add(str + "$");
      m_leaves = new ArrayList<Node>();
   }

   public void addString(String str) {
      m_sequences.add(str + "$");
   }

   public String substr(Position start, Position end) {
      if (start == null)
          System.out.println("start is null!");
      if (end == null)
          System.out.println("end is null!");

      System.out.println(Integer.toString(start.idx) + " : " + Integer.toString(start.position) + " : " + 
                         Integer.toString(end.position) + " size=" + Integer.toString(m_sequences.get(start.idx).length()));
      return m_sequences.get(start.idx).substring(start.position, end.position + 1); // end is inclusive
   }

   public Node buildTree() {
      if (root == null)
         root = new Node(Node.TYPE_ROOT);
      for (int idx = 0; idx < m_sequences.size(); idx++) {
         for(int i = 0; i < m_sequences.get(idx).length(); i++) {
            insert(root,i, i, m_sequences.get(idx).length() - 1, idx);
         }
      }
      return root;
   }

   public void insert(Node n, int pos, int start, int end, int idx) {
      int count = 0;
      //System.out.println("before while");
      for(int i = 0; i < n.getNumChildren(); i++) {
          //System.out.println("i : " + i + " " + m_sequences.get(idx).charAt(start) + " ? " + m_sequences.get(n.getChildAt(i).getStart().idx).charAt(n.getChildAt(i).getStart().position));
          if (m_sequences.get(n.getChildAt(i).getStart().idx).charAt(n.getChildAt(i).getStart().position) != '$' &&
              (m_sequences.get(n.getChildAt(i).getStart().idx).charAt(n.getChildAt(i).getStart().position) == m_sequences.get(idx).charAt(start) ||
              (n.getChildAt(i).getStart().idx != idx && m_sequences.get(n.getChildAt(i).getStart().idx).charAt(n.getChildAt(i).getStart().position) == 'N') ||
              m_sequences.get(idx).charAt(start) == 'N')) {
             count++;
             int j = 1;
             while(j <= n.getChildAt(i).getEnd().position - n.getChildAt(i).getStart().position && j <= end - start
             && (m_sequences.get(n.getChildAt(i).getStart().idx).charAt(n.getChildAt(i).getStart().position + j) == m_sequences.get(idx).charAt(start + j)
             || (n.getChildAt(i).getStart().idx != idx && m_sequences.get(n.getChildAt(i).getStart().idx).charAt(n.getChildAt(i).getStart().position + j) == 'N')
             || m_sequences.get(idx).charAt(start + j) == 'N')) {
                j++;
             } //compare edge label

             j--;
             if (j == n.getChildAt(i).getEnd().position - n.getChildAt(i).getStart().position) { //entire label match
                if (j + 1 < end - start && n.getChildAt(i).getNumChildren() > 0) {//Node has children
                   //System.out.println("recursive insert!");
                   insert(n.getChildAt(i), pos, start + j + 1, end, idx);
                }
                else {
                   /*System.out.println(substr(new Position(idx,start + j + 1), 
                                                                     new Position(idx,end)) + " //leafnode");*/
                   Node inode = new Node(new Position(idx,start), new Position(idx,end-1), Node.TYPE_INTERNAL,
                                                                     new Position(idx,pos)); //add new internal node
                   /*Node inode = new Node(n.getChildAt(i).getStart(), new Position(n.getChildAt(i).getStart().idx,n.getChildAt(i).getStart().position + j - 1),
                                                                     Node.TYPE_INTERNAL, new Position(idx,pos)); //add new internal node*/
                   Node node = new Node(new Position(idx,start + j), new Position(idx,end), Node.TYPE_LEAF, 
                                                                     new Position(idx,pos)); //add new leaf with $
                   node.setParent(inode);
                   inode.addChild(node);
                   inode.addPosition(n.getChildAt(i).getPositionAt(0));
                   inode.setParent(n);
                   n.getChildAt(i).setParent(inode);
                   n.getChildAt(i).resetStart(j);
                   inode.addChild(n.getChildAt(i));
                   Node rnode = n.removeChildAt(i);
                   //System.out.println("removed1 : " + substr(rnode.getStart(),rnode.getEnd()));
                   propagatePositionUp(n, new Position(idx,pos));
                   n.addChild(i, inode);
                   m_leaves.add(node);
                }
            }
            else {
               //System.out.println(start + " //inode : " + j);
               Node inode = new Node(new Position(idx,start), new Position(idx,start + j), Node.TYPE_INTERNAL); //new internal n
               /*Node inode = new Node(new Position(n.getChildAt(i).getStart().idx,n.getChildAt(i).getStart().position), 
                            new Position(n.getChildAt(i).getStart().idx,n.getChildAt(i).getStart().position + j), Node.TYPE_INTERNAL);*/
               Node lnode = new Node(new Position(idx,start + j + 1), new Position(idx,end), Node.TYPE_LEAF, new Position(idx,pos)); //create new leaf lnode
               lnode.setParent(inode);
               m_leaves.add(lnode); //add to list of leaves
               inode.addChild(lnode);
               inode.setParent(n);
               n.getChildAt(i).resetStart(j + 1);
               n.getChildAt(i).setParent(inode);
               inode.addChild(n.getChildAt(i));
               inode.addPosition(new Position(idx,pos));
               for(Position p : n.getChildAt(i).getPositions()) {
                  inode.addPosition(p);
               }
               Node rnode = n.removeChildAt(i);
               //System.out.println("removed2 : " + substr(rnode.getStart(),rnode.getEnd()));
               n.addChild(i, inode);
               propagatePositionUp(n, new Position(idx,pos));
            }
         }
      }
      if (count == 0) {
         //System.out.println("count 0 : " + m_sequences.get(idx).charAt(start));
         //System.out.println(substr(new Position(idx,start), new Position(idx,end)) + " //leafnode");
         Node node = new Node(new Position(idx,start), new Position(idx,end), Node.TYPE_LEAF, new Position(idx,pos));
         node.setParent(n);
         n.addChild(node);
         //n.addPosition(new Position(idx,pos));
         propagatePositionUp(n, new Position(idx,pos));
         m_leaves.add(node);  //add to list of leaves
      }

   }

   // add position to all nodes upward in the tree
   public void propagatePositionUp(Node n, Position p) {
      Node currNode = n;
      do {
         currNode.addPosition(p);
         currNode = currNode.getParent();
      } while(currNode != null);
   }

   public ArrayList<RepeatSeq> findRepeats(int len, double multiply) { //bottom up search
      ArrayList<RepeatSeq> repeats = new ArrayList<RepeatSeq>();
      ArrayList<Long> visitedNodes = new ArrayList<Long>();
      String label;
      Node currNode;
      double expected;
      boolean isShort;
      for (Node leaf : m_leaves) {
         label = "";
         isShort = false;
         currNode = leaf.getParent();
         while (!isShort && currNode.getType() == Node.TYPE_INTERNAL && !isVisited(visitedNodes, currNode)) {
            label = getLabel(currNode);
            if (label.length() >= len) {
               expected = calcExpected(label);
               //System.out.println(Float.toString(expected));
               if (multiply > 0.0){
                  if (currNode.getNumPositions() >= expected * multiply) {
                     //System.out.println("findRepeats: " + Float.toString(expected));
                     repeats.add(new RepeatSeq(label, expected, currNode.getNumPositions(),currNode.getPositions()));
                  }
               }
               else {
                  repeats.add(new RepeatSeq(label, expected, currNode.getNumPositions(),currNode.getPositions()));
               }
            }
            else {
               isShort = true; // label is too short, so get out of the loop.
            }
            visitedNodes.add(currNode.getId());
            currNode = currNode.getParent();
         }
      }
      return repeats;
   }

   public String getLabel(Node node) {
      int i = 0;
      int idx = node.getPositionAt(i).idx;
      while(i < node.getNumPositions() && idx != node.getStart().idx) {
         i++;
         idx = node.getPositionAt(i).idx;
      }
      return substr(node.getPositionAt(i), node.getEnd());
   }

   public boolean isVisited(ArrayList<Long> list, Node node) {
      long id = node.getId();
      for(long l : list) {
         if(l == id)
            return true;
      }
      return false;
   }

   public void calcNucleotideContent() {
      totalA = 0.0f;
      totalC = 0.0f;
      totalG = 0.0f;
      totalT = 0.0f;
      totalNucleotides = 0;

      for (String str : m_sequences) {
         for (char c : str.toCharArray()) {
            switch(c) {
            case 'A' :
               totalA += 1.0;
               totalNucleotides++;
               break;
            case 'C' :
               totalC += 1.0;
               totalNucleotides++;
               break;
            case 'G' :
               totalG += 1.0;
               totalNucleotides++;
               break;
            case 'T' :
               totalT += 1.0;
               totalNucleotides++;
               break;
            case 'N' :
               totalA += 0.25;
               totalC += 0.25;
               totalG += 0.25;
               totalT += 0.25;
               totalNucleotides++;
               break;
            }
         }
      }

   }

   public double calcExpected(String str) {
      double a, t, c, g;
      double aPer = totalA / (double)totalNucleotides;
      double tPer = totalT / (double)totalNucleotides;
      double cPer = totalC / (double)totalNucleotides;
      double gPer = totalG / (double)totalNucleotides;
      a = t = c = g = 0.0;
      for(char l : str.toCharArray()) {
         switch(l){
         case 'A': a++; break;
         case 'T': t++; break;
         case 'C': c++; break;
         case 'G': g++; break;
         }
      }
      //System.out.println(Float.toString(pow(aPer,a) * pow(tPer,t) * pow(cPer,c) * pow(gPer,g) * totalNucleotides);
      
      return Math.pow(aPer,a) * Math.pow(tPer,t) * Math.pow(cPer,c) * Math.pow(gPer,g) * totalNucleotides;
   }

   public static void main(String[] args){
      String str;
      try{
         str = ReadGFF.readFasta("bioGUI/ananase_contig.txt");
         str = "TTTACTTTACGG";
         SuffixTree sf = new SuffixTree(str);
         sf.addString("ANGG");
         Node r = sf.buildTree();
         System.out.println(r.getNumChildren());
         for (int i =0; i < r.getNumChildren(); i++) {
            System.out.println(i + " : " + sf.substr(r.getChildAt(i).getStart(),r.getChildAt(i).getEnd()));
            System.out.println("length : " + r.getChildAt(i).getLength());
            for(Position p : r.getChildAt(i).getPositions())
               System.out.println("position : " + p.toString());
            for (Node n : r.getChildAt(i).getChildren()) {
               System.out.println(sf.substr(n.getStart(),n.getEnd()));
               System.out.println("length : " + (r.getChildAt(i).getLength() + n.getLength()));
               for(Position p : n.getPositions())
                  System.out.println("       position : " + p.toString());
            }
         }
         System.out.println("number of leaves = " + Integer.toString(sf.m_leaves.size()));
         sf.calcNucleotideContent();
         ArrayList<RepeatSeq> repeats =  sf.findRepeats(2, 0.0f);
         System.out.println(Integer.toString(repeats.size()));
         for (RepeatSeq rep : repeats) {
            System.out.println(rep.toString());
         }
      } catch (IOException e) {
         System.err.println(e.getMessage());
      }

   }
}

import CCfg.*;

public class GetConfigParser {
	
	public void parse(String fileName)
	{
	      DocumentFacade mydoc = DocumentFacade.loadDoc(fileName);
	      CConfigFacade conf = mydoc.getConfig();

	      int shapingTime_ns = conf.getInt("Node[CoBo].Instance[0].AsAd[0].Aget[0].Global.Reg1.peackingTime");
	      System.out.println("Shaping time: " + shapingTime_ns + " ns");
	      
	      boolean isSelected = conf.getBool("Node[CoBo].Instance[0].AsAd[0].Aget[0].channel[3].isSelectedforTestMode");
	      System.out.println("isSelected: " + isSelected);
	      
	      System.out.println();
	      
	      StringVector nodes = new StringVector();
	      conf.enumNode("Node[CoBo]", nodes);
	      for (int i=0; i < nodes.size(); ++i)
	      {
	    	  System.out.println(nodes.get(i));
	      }
	      
	      System.out.println();
	      
	      nodes.clear();
	      conf.enumNode("Node[CoBo].Instance[0]", nodes);
	      for (int i=0; i < nodes.size(); ++i)
	      {
	    	  if (! nodes.get(i).startsWith("AsAd")) continue;
	    	  System.out.println(nodes.get(i));
	      }
	}
	
	public static void main(String[] args)
	{
		if (args.length <= 0)
			System.exit(1);
		GetConfigParser parser = new GetConfigParser();
		parser.parse(args[0]);
	}
}
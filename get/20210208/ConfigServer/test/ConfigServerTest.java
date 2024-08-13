import org.junit.* ;
import static org.junit.Assert.* ;
import CCfg.*;

public class ConfigServerTest {
   @Test
   public void test_loadDoc() {
      System.out.println("loadDoc() Test") ;
      DocumentFacade myDoc = DocumentFacade.loadDoc("ConfigTest.xml");  
	  CConfigFacade conf = myDoc.getConfig();
	  assertEquals(conf.getName(""), "Test");
   }
   
	@Test
   public void test_newDoc() {
      System.out.println("newDoc() Test") ;
      DocumentFacade myDoc = DocumentFacade.loadDoc("ConfigTest.xml");  
	  CConfigFacade conf = myDoc.getConfig();
	  assertEquals(conf.getName(""), "Test");
   }
}

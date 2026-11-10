using System.Windows;

namespace Cidep.App;

public partial class App : Application
{
    protected override void OnStartup(StartupEventArgs args)
    {
        base.OnStartup(args);
        var window = new MainWindow();
        MainWindow = window;
        if (args.Args.Contains("--smoke-test"))
            window.Loaded += async (_, _) =>
            {
                try
                {
                    await window.SmokeTest();
                    Shutdown(0);
                }
                catch (Exception error)
                {
                    System.IO.Directory.CreateDirectory("artifacts");
                    System.IO.File.WriteAllText("artifacts/wpf-error.txt", error.ToString());
                    Shutdown(1);
                }
            };
        window.Show();
    }
}
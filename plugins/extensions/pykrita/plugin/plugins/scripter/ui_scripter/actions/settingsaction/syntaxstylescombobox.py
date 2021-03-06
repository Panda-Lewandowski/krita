from PyQt5.QtWidgets import QComboBox
from scripter.ui_scripter.syntax import syntaxstyles


class SyntaxStylesComboBox(QComboBox):

    def __init__(self, highlight, parent=None):
        super(SyntaxStylesComboBox, self).__init__(parent)

        self.highlight = highlight
        self.styleClasses = [syntaxstyles.DefaultSyntaxStyle, syntaxstyles.PythonVimSyntaxStyle]

        for styleClass in self.styleClasses:
            className = styleClass.__name__
            self.addItem(className)

            if className == type(self.highlight.getSyntaxStyle()).__name__:
                self.setCurrentIndex(self.findText(className))

        self.currentIndexChanged.connect(self._currentIndexChanged)

    def _currentIndexChanged(self, index):
        self.highlight.setSyntaxStyle(getattr(syntaxstyles, self.itemText(index))())
        self.highlight.rehighlight()

    def readSettings(self, settings):
        syntaxStyle = settings.value('syntaxStyle', '')

        if syntaxStyle:
            self.setCurrentIndex(self.findText(syntaxStyle))

    def writeSettings(self, settings):
        settings.setValue('syntaxStyle', type(self.highlight.getSyntaxStyle()).__name__)
